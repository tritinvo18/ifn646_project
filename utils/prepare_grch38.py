#!/usr/bin/env python3
"""Download and prepare NCBI GCF_000001405.26 GRCh38 reference FASTA.

Creates under <project-root>/data by default:
  GCF_000001405.26_GRCh38_genomic.fna.gz
  GRCh38.fa
  GRCh38.fa.fai

1. Run it: 

python utils/prepare_grch38.py \
    --project-root .

2. Delete archive after FASTA and index have been validated:

python scripts/prepare_grch38_reference.py \
    --project-root . \
    --delete-archive

3. Replace an existing incorrect FASTA

python scripts/prepare_grch38_reference.py \
  --project-root . \
  --force-decompress \
  --force-index

4. Force a completely new download

python scripts/prepare_grch38_reference.py \
  --project-root . \
  --force-download \
  --force-decompress \
  --force-index

"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ACCESSION = "GCF_000001405.26"
ASSEMBLY = "GRCh38"
ARCHIVE_NAME = f"{ACCESSION}_{ASSEMBLY}_genomic.fna.gz"
BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
    f"{ACCESSION}_{ASSEMBLY}"
)
DOWNLOAD_URL = f"{BASE_URL}/{ARCHIVE_NAME}"
MD5_URL = f"{BASE_URL}/md5checksums.txt"

EXPECTED_PRIMARY = {
    "NC_000001.11",
    "NC_000002.12",
    "NC_000003.12",
    "NC_000004.12",
    "NC_000005.10",
    "NC_000006.12",
    "NC_000007.14",
    "NC_000008.11",
    "NC_000009.12",
    "NC_000010.11",
    "NC_000011.10",
    "NC_000012.12",
    "NC_000013.11",
    "NC_000014.9",
    "NC_000015.10",
    "NC_000016.10",
    "NC_000017.11",
    "NC_000018.10",
    "NC_000019.10",
    "NC_000020.11",
    "NC_000021.9",
    "NC_000022.11",
    "NC_000023.11",
    "NC_000024.10",
}


def human_bytes(value: int) -> str:
    size = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{value} B"


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "task3-reference-preparer/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8")


def expected_md5() -> str:
    text = fetch_text(MD5_URL)
    pattern = re.compile(rf"^([0-9a-fA-F]{{32}})\s+\./{re.escape(ARCHIVE_NAME)}$", re.MULTILINE)
    match = pattern.search(text)
    if not match:
        raise RuntimeError(f"Could not locate {ARCHIVE_NAME} in {MD5_URL}")
    return match.group(1).lower()


def file_md5(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


def download_with_resume(url: str, destination: Path) -> None:
    partial = Path(str(destination) + ".part")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": "task3-reference-preparer/1.0"}
    if existing:
        headers["Range"] = f"bytes={existing}-"

    request = urllib.request.Request(url, headers=headers)
    try:
        response = urllib.request.urlopen(request, timeout=120)
    except urllib.error.HTTPError as exc:
        if exc.code == 416 and partial.exists():
            partial.replace(destination)
            return
        raise

    status = getattr(response, "status", response.getcode())
    if existing and status != 206:
        print("Server did not honour resume request; restarting download.")
        existing = 0
        partial.unlink(missing_ok=True)

    content_length = response.headers.get("Content-Length")
    total = existing + int(content_length) if content_length else None
    mode = "ab" if existing else "wb"
    downloaded = existing
    last_report = -1

    print(f"Downloading: {url}")
    if existing:
        print(f"Resuming at {human_bytes(existing)}")

    with response, partial.open(mode) as output:
        while True:
            chunk = response.read(8 * 1024 * 1024)
            if not chunk:
                break
            output.write(chunk)
            downloaded += len(chunk)
            if total:
                percent = int(downloaded * 100 / total)
                if percent >= last_report + 2:
                    print(
                        f"  {percent:3d}%  {human_bytes(downloaded)} / {human_bytes(total)}",
                        flush=True,
                    )
                    last_report = percent
            elif downloaded // (128 * 1024 * 1024) > last_report:
                print(f"  downloaded {human_bytes(downloaded)}", flush=True)
                last_report = downloaded // (128 * 1024 * 1024)

    partial.replace(destination)
    print(f"Downloaded {human_bytes(destination.stat().st_size)} to {destination}")


def verify_gzip(path: Path) -> None:
    with gzip.open(path, "rb") as handle:
        first = handle.read(1)
    if first != b">":
        raise RuntimeError(f"Archive does not decompress to FASTA: {path}")


def decompress_atomic(archive: Path, fasta: Path, force: bool) -> None:
    if fasta.exists() and not force:
        print(f"FASTA exists; retaining: {fasta}")
        return

    temporary = Path(str(fasta) + ".tmp")
    temporary.unlink(missing_ok=True)
    print(f"Decompressing {archive.name} to {fasta}")
    try:
        with gzip.open(archive, "rb") as source, temporary.open("wb") as target:
            shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
        if temporary.stat().st_size == 0:
            raise RuntimeError("Decompressed FASTA is empty")
        temporary.replace(fasta)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    print(f"Created {fasta} ({human_bytes(fasta.stat().st_size)})")


def fasta_contigs(path: Path) -> set[str]:
    contigs: set[str] = set()
    with path.open("rt", encoding="ascii") as handle:
        first_line = handle.readline()
        if not first_line.startswith(">"):
            raise RuntimeError(f"Not a FASTA file: {path}")
        contigs.add(first_line[1:].split()[0])
        for line in handle:
            if line.startswith(">"):
                contigs.add(line[1:].split()[0])
    return contigs


def validate_fasta(path: Path) -> None:
    contigs = fasta_contigs(path)
    missing = sorted(EXPECTED_PRIMARY - contigs)
    if missing:
        raise RuntimeError(f"FASTA is missing expected primary chromosomes: {missing}")
    print(f"Validated FASTA with {len(contigs)} sequence records.")
    print("Primary chromosome accessions include NC_000001.11 through NC_000024.10.")


def build_index(fasta: Path, force: bool) -> Path:
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools is required but was not found in PATH")
    index = Path(str(fasta) + ".fai")
    if force:
        index.unlink(missing_ok=True)
    if index.exists() and index.stat().st_mtime >= fasta.stat().st_mtime:
        print(f"FASTA index is current: {index}")
        return index
    print(f"Building FASTA index: {index}")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    if not index.exists() or index.stat().st_size == 0:
        raise RuntimeError(f"samtools did not create a valid index: {index}")
    return index


def validate_index(index: Path) -> None:
    names = set()
    with index.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise RuntimeError(f"Malformed FASTA index row: {line.rstrip()}")
            names.add(fields[0])
    missing = sorted(EXPECTED_PRIMARY - names)
    if missing:
        raise RuntimeError(f"FASTA index is missing primary chromosomes: {missing}")
    print(f"Validated index with {len(names)} sequence records: {index}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--force-decompress", action="store_true")
    parser.add_argument("--force-index", action="store_true")
    parser.add_argument(
        "--delete-archive",
        action="store_true",
        help="Delete the .fna.gz only after FASTA and index validation succeeds",
    )
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)

    archive = data / ARCHIVE_NAME
    fasta = data / "GRCh38.fa"

    if args.force_download:
        archive.unlink(missing_ok=True)
        Path(str(archive) + ".part").unlink(missing_ok=True)

    checksum = expected_md5()
    print(f"Expected NCBI MD5: {checksum}")

    if archive.exists():
        observed = file_md5(archive)
        if observed == checksum:
            print(f"Archive already exists and checksum matches: {archive}")
        else:
            print("Existing archive checksum is incorrect; downloading again.")
            archive.unlink()
            download_with_resume(DOWNLOAD_URL, archive)
    else:
        download_with_resume(DOWNLOAD_URL, archive)

    observed = file_md5(archive)
    if observed != checksum:
        raise RuntimeError(
            f"MD5 mismatch for {archive}: expected {checksum}, observed {observed}"
        )
    print("Archive MD5 checksum passed.")

    verify_gzip(archive)
    decompress_atomic(archive, fasta, args.force_decompress)
    validate_fasta(fasta)
    index = build_index(fasta, args.force_index)
    validate_index(index)

    if args.delete_archive:
        archive.unlink()
        print(f"Deleted compressed archive: {archive}")

    print("\nReference preparation complete:")
    print(f"  FASTA: {fasta}")
    print(f"  index: {index}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted. The partial download is retained for resuming.", file=sys.stderr)
        raise SystemExit(130)
