import re
import gzip
import pandas as pd
from pathlib import Path

def attrs(text):
    out = {}
    for item in text.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1); out[k] = v
    return out

def transcript_match(a, accession):
    values = "|".join(map(str, a.values()))
    return accession in values or accession.split(".")[0] in values

def collect_exons(path, chrom_accession, preferred_transcript, gene=None):
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip().split("\t")
            if len(p) != 9 or p[0] != chrom_accession or p[2] != "exon": continue
            a = attrs(p[8])
            if gene:
                if gene not in {a.get("gene"), a.get("Name"), a.get("gene_name")} and f"gene={gene}" not in p[8]: continue
            if not transcript_match(a, preferred_transcript): continue
            rows.append({"seqid":p[0], "start":int(p[3]), "end":int(p[4]), "strand":p[6]})
    return pd.DataFrame(rows).drop_duplicates()

def uncompress_fasta(path):
    path = Path(path)
    if path.suffix != ".gz": return path
    out = path.with_suffix("")
    if not out.exists():
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            while block := src.read(1024*1024): dst.write(block)
    return out

def rc(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def load_crispor(path):
    engine = "xlrd" if str(path).lower().endswith(".xls") else "openpyxl"
    preview = pd.read_excel(path, sheet_name="guides", header=None, engine=engine)
    rows = preview.index[preview.iloc[:, 0].astype(str).str.strip().eq("#guideId")].tolist()
    if len(rows) != 1:
        raise ValueError(f"Expected one #guideId header; found {rows}")
    d = pd.read_excel(path, sheet_name="guides", header=rows[0], engine=engine)
    d.columns = d.columns.astype(str).str.strip()
    for old in ["Doench '16-Score", "Doench 2016-Score", "Doench2016"]:
        if old in d.columns:
            d = d.rename(columns={old: "Efficiency"})
            break
    d["targetSeq"] = d["targetSeq"].astype(str).str.upper().str.replace(r"[^ACGTN]", "", regex=True)
    for c in ["mitSpecScore", "cfdSpecScore", "Efficiency"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["#guideId", "targetSeq"]).copy()

def all_starts(sequence, motif):
    return [m.start() for m in re.finditer(f"(?={motif})", sequence)]

def map_guide(row, search_sequence, first_exon):
    target = row["targetSeq"]
    if len(target) < 23:
        return pd.Series({
            "mapping_status": "INVALID_LENGTH", "mapped_orientation": pd.NA,
            "protospacer": target[:20], "protospacer_start0": pd.NA,
            "protospacer_end0": pd.NA, "cut_boundary0": pd.NA,
            "overlaps_first_exon": False, "cut_inside_first_exon": False,
        })
    target23 = target[:23]
    candidates = []
    for orientation, genomic23 in [("forward", target23), ("reverse", rc(target23))]:
        for start0 in all_starts(search_sequence, genomic23):
            if orientation == "forward":
                protospacer_start0 = start0
                protospacer_end0 = start0 + 20
                cut_boundary0 = start0 + 17
                protospacer = target23[:20]
            else:
                protospacer_start0 = start0 + 3
                protospacer_end0 = start0 + 23
                cut_boundary0 = start0 + 6
                protospacer = target23[:20]
            candidates.append({
                "mapped_orientation": orientation,
                "protospacer": protospacer,
                "protospacer_start0": protospacer_start0,
                "protospacer_end0": protospacer_end0,
                "cut_boundary0": cut_boundary0,
            })
    if len(candidates) != 1:
        return pd.Series({
            "mapping_status": "NOT_FOUND" if not candidates else "AMBIGUOUS",
            "mapped_orientation": pd.NA, "protospacer": target23[:20],
            "protospacer_start0": pd.NA, "protospacer_end0": pd.NA,
            "cut_boundary0": pd.NA, "overlaps_first_exon": False,
            "cut_inside_first_exon": False,
        })
    hit = candidates[0]
    
    # Check overlaps with first_exon inside the search_sequence coordinate space
    search_start = max(1, int(first_exon["start"]) - 25)
    exon_start0 = int(first_exon["start"]) - search_start
    exon_end0 = int(first_exon["end"]) - search_start + 1
    
    overlap = max(0, min(hit["protospacer_end0"], exon_end0) - max(hit["protospacer_start0"], exon_start0))
    cut_inside = (exon_start0 <= hit["cut_boundary0"] < exon_end0)
    
    return pd.Series({
        "mapping_status": "UNIQUE", "mapped_orientation": hit["mapped_orientation"],
        "protospacer": hit["protospacer"], "protospacer_start0": hit["protospacer_start0"],
        "protospacer_end0": hit["protospacer_end0"], "cut_boundary0": hit["cut_boundary0"],
        "overlaps_first_exon": overlap > 0,
        "cut_inside_first_exon": cut_inside
    })
