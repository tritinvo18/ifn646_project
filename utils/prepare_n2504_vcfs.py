#!/usr/bin/env python3
"""
Prepare compact Task 3 regional VCFs for N=2,504 cohort by remote indexed 1000G queries.

Run from project root:
  python utils/prepare_n2504_vcfs.py --chromosomes target
  or
  python utils/prepare_n2504_vcfs.py --chromosomes autosomes

Outputs are saved in:
  data/vcfs/n2504/sampled_chr{chrom}.norm.vcf.gz
"""
from __future__ import annotations
import argparse, csv, os, re, shutil, subprocess, sys
from pathlib import Path

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
TEMPLATE = BASE + "/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

REFSEQ = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
    "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
    "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
    "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11"
}

ORDER = {str(i): i for i in range(1, 23)} | {"X": 23}

ON_TARGET = {
    "4": [(76035800 - 1, 76036350, "CXCL11_on_target")],
    "18": [(63887550 - 1, 63887950, "SERPINB2_on_target")]
}

def get_bcftools_bin() -> str:
    path = shutil.which("bcftools")
    if path:
        return path
    candidate = "/home/jimmy/.conda/envs/ifn646/bin/bcftools"
    if os.path.exists(candidate):
        return candidate
    raise SystemExit("ERROR: bcftools binary not found in PATH or conda env")

BCFTOOLS = get_bcftools_bin()

def run(cmd, capture=False):
    print("+", " ".join(map(str, cmd)), flush=True)
    r = subprocess.run(list(map(str, cmd)), check=True, text=True, stdout=subprocess.PIPE if capture else None)
    return r.stdout if capture else ""

def canon(x):
    x = x.strip()
    m = re.fullmatch(r"chr([0-9]{1,2}|X|Y)", x)
    if m:
        x = m.group(1)
    if x in REFSEQ:
        return x
    for c, a in REFSEQ.items():
        if x == a:
            return c
    return None

def samples(path):
    s = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    if len(s) != 2504 or len(set(s)) != 2504:
        raise ValueError(f"Expected 2,504 unique sample IDs in {path}, found {len(set(s))}")
    return s

def regions(path):
    out = {}
    excluded = set()
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        c = canon(f[0])
        if c is None:
            excluded.add(f[0])
            continue
        a, b = int(f[1]), int(f[2])
        label = f[3] if len(f) > 3 else f"site_{n}"
        if a < 0 or b <= a:
            raise ValueError(f"Invalid BED interval line {n}")
        out.setdefault(c, []).append((a, b, label))
    for c, v in ON_TARGET.items():
        out.setdefault(c, []).extend(v)
    return {c: sorted(set(v)) for c, v in out.items()}, excluded

def remote_contig(url, chrom):
    h = run([BCFTOOLS, "view", "-h", url], True)
    contigs = set(re.findall(r"##contig=<ID=([^,>]+)", h))
    for x in [chrom, f"chr{chrom}", REFSEQ[chrom]]:
        if x in contigs:
            return x
    raise ValueError(f"Cannot resolve chromosome {chrom}; examples={sorted(contigs)[:10]}")

def ready(path, expected, refseq):
    if not path.exists() or not any(Path(str(path) + x).exists() for x in [".csi", ".tbi"]):
        return False
    try:
        observed = set(run([BCFTOOLS, "query", "-l", path], True).splitlines())
        header = run([BCFTOOLS, "view", "-h", path], True)
        return observed == set(expected) and re.search(rf"##contig=<ID={re.escape(refseq)}(?:,|>)", header) is not None
    except subprocess.CalledProcessError:
        return False

def delete_files(paths):
    count = 0
    for p in paths:
        p = Path(p)
        if p.is_file():
            p.unlink()
            count += 1
    return count

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--chromosomes", choices=["target", "autosomes", "all"], default="target",
                    help="'target' extracts Chr 4 & 18 (on-target genes & off-targets), 'autosomes' extracts Chr 1-22")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--keep-temporary", action="store_true")
    a = ap.parse_args()

    root = Path(a.project_root).resolve()
    data = root / "data"
    tables = root / "results" / "tables"
    work = data / "offtarget_regions_n2504"
    work.mkdir(parents=True, exist_ok=True)

    ref = data / "GRCh38.fa"
    sample_file = data / "n2504" / "sample_ids.txt"
    bed_all = tables / "task3_offtarget_regions.bed"

    for p in [ref, Path(str(ref) + ".fai"), sample_file, bed_all]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required file: {p}")

    expected = samples(sample_file)
    reg, excluded = regions(bed_all)

    if a.chromosomes == "target":
        wanted = {"4", "18"}
    elif a.chromosomes == "autosomes":
        wanted = {str(i) for i in range(1, 23)}
    else:
        wanted = {str(i) for i in range(1, 23)} | {"X"}

    selected = sorted(set(reg) & wanted, key=lambda x: ORDER[x])
    print(f"Selected chromosomes for N=2504 extraction: {','.join(selected)}")
    if excluded:
        print("Excluded non-primary contigs:", ",".join(sorted(excluded)))

    vcf_dir = data / "vcfs" / "n2504"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    manifest = []
    for chrom in selected:
        print(f"\n=== Chromosome {chrom} (N=2,504) ===")
        output = vcf_dir / f"sampled_chr{chrom}.norm.vcf.gz"
        marker = Path(str(output) + ".n2504_regions_complete")
        url = TEMPLATE.format(chrom=chrom)

        if not a.force and marker.exists() and ready(output, expected, REFSEQ[chrom]):
            print("SKIP: complete output exists")
            manifest.append([chrom, "SKIPPED", str(output), url, ""])
            continue

        raw = work / f"sampled_chr{chrom}.raw.bcf"
        renamed = work / f"sampled_chr{chrom}.renamed.bcf"
        bed = work / f"task3_chr{chrom}.bed"
        rename = work / f"rename_chr{chrom}.tsv"
        tmp = Path(str(output) + ".tmp.vcf.gz")

        try:
            sc = remote_contig(url, chrom)
            print(f"Remote contig: {sc}")
            bed.write_text("".join(f"{sc}\t{x}\t{y}\t{z}\n" for x, y, z in reg[chrom]))
            rename.write_text(f"{sc}\t{REFSEQ[chrom]}\n")

            run([BCFTOOLS, "view", "-S", sample_file, "-R", bed, "-Ob", "-o", raw, url])
            run([BCFTOOLS, "index", "-f", raw])
            run([BCFTOOLS, "annotate", "--rename-chrs", rename, "-Ob", "-o", renamed, raw])
            run([BCFTOOLS, "index", "-f", renamed])
            run([BCFTOOLS, "norm", "-f", ref, "-c", "e", "-m", "-any", "-Oz", "-o", tmp, renamed])
            run([BCFTOOLS, "index", "-f", tmp])

            delete_files([output, Path(str(output) + ".csi"), Path(str(output) + ".tbi")])
            tmp.rename(output)
            for suf in [".csi", ".tbi"]:
                q = Path(str(tmp) + suf)
                if q.exists():
                    q.rename(Path(str(output) + suf))

            if not ready(output, expected, REFSEQ[chrom]):
                raise RuntimeError(f"Output verification failed: {output}")

            n = int(run([BCFTOOLS, "index", "--nrecords", output], True).strip() or 0)
            print(f"Retained {n} regional variant records across 2,504 individuals.")
            marker.write_text(f"source={url}\nregions={len(reg[chrom])}\nrecords={n}\nsamples=2504\n")

            if not a.keep_temporary:
                delete_files([bed, rename, raw, renamed, Path(str(raw) + ".csi"), Path(str(renamed) + ".csi")])
            manifest.append([chrom, "CREATED", str(output), url, ""])
        except Exception as e:
            print(f"ERROR chromosome {chrom}: {e}", file=sys.stderr)
            manifest.append([chrom, "FAILED", str(output), url, str(e)])

    manifest_file = tables / "n2504" / "task3_n2504_vcf_manifest.tsv"
    manifest_file.parent.mkdir(parents=True, exist_ok=True)
    with manifest_file.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["chromosome", "status", "output_vcf", "source", "error"])
        w.writerows(manifest)

    if not a.keep_temporary:
        delete_files(work.glob("*"))
        try:
            work.rmdir()
        except OSError:
            pass

    failed = [x for x in manifest if x[1] == "FAILED"]
    print(f"\nManifest saved to: {manifest_file}")
    if failed:
        raise SystemExit(f"Completed with {len(failed)} failed chromosome(s).")
    print("N=2,504 Regional VCF preparation complete!")

if __name__ == "__main__":
    main()
