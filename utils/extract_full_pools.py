#!/usr/bin/env python3
"""
Extract the full pool of eligible CRISPR guide RNAs for each chromosome separately
from Task 2 audit logs, without modifying the original Task 2 notebooks.

Chromosome 4: CXCL11 (CRISPRi) - all guides passing hard eligibility filters
Chromosome 18: SERPINB2 (CRISPRa) - all guides passing hard eligibility filters
"""

from pathlib import Path
import pandas as pd

TABLE_DIR = Path(__file__).resolve().parent.parent / "results" / "tables"

def extract_chromosome_full_pool(gene: str, mode: str, chrom: str) -> pd.DataFrame:
    audit_file = TABLE_DIR / f"task2_{gene.lower()}_{mode.lower()}_all_guides_audit.csv"
    if not audit_file.exists():
        raise FileNotFoundError(f"Audit file not found: {audit_file}")
    
    df = pd.read_csv(audit_file)
    eligible_df = df.loc[df["eligible"]].copy()
    eligible_df["gene"] = gene
    eligible_df["mode"] = mode
    eligible_df["chrom"] = chrom
    
    out_file = TABLE_DIR / f"task2_{gene.lower()}_{mode.lower()}_full_pool.csv"
    eligible_df.to_csv(out_file, index=False)
    print(f"[{gene} - Chr {chrom}] Extracted {len(eligible_df)} eligible guides -> {out_file.name}")
    return eligible_df

if __name__ == "__main__":
    print("Extracting full pools for each chromosome separately...")
    cxcl11_full = extract_chromosome_full_pool("CXCL11", "CRISPRi", "4")
    serpinb2_full = extract_chromosome_full_pool("SERPINB2", "CRISPRa", "18")
    print(f"Total eligible guides across both chromosomes: {len(cxcl11_full) + len(serpinb2_full)}")
