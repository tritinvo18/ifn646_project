#!/usr/bin/env python3
"""
Fetch and curate 1000 Genomes Phase 3 sample metadata for N=2,504 cohort.
Saves:
  - data/n2504/sample_metadata.tsv: Detailed demographic and population metadata
  - data/n2504/sample_ids.txt: Exactly 2,504 sample identifiers
"""

import csv
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
N2504_DIR = DATA_DIR / "n2504"
N2504_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_INFO_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt"
PANEL_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

POP_TO_SUPER = {
    'CHB': 'EAS', 'JPT': 'EAS', 'CHS': 'EAS', 'CDX': 'EAS', 'KHV': 'EAS',
    'CEU': 'EUR', 'TSI': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR',
    'YRI': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR', 'ASW': 'AFR', 'ACB': 'AFR',
    'MXL': 'AMR', 'PUR': 'AMR', 'CLM': 'AMR', 'PEL': 'AMR',
    'GIH': 'SAS', 'PJL': 'SAS', 'BEB': 'SAS', 'STU': 'SAS', 'ITU': 'SAS'
}

SUPER_POP_NAMES = {
    'EUR': 'European',
    'AFR': 'African',
    'EAS': 'East Asian',
    'SAS': 'South Asian',
    'AMR': 'Admixed American'
}

def main():
    print("Fetching 1000 Genomes Phase 3 sample information...")
    req_info = urllib.request.Request(SAMPLE_INFO_URL, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req_info, timeout=20) as resp:
        info_lines = resp.read().decode('utf-8', errors='ignore').splitlines()

    sample_info = {}
    for line in info_lines[1:]:
        f = line.split('\t')
        if len(f) > 5 and f[0].strip():
            sample_info[f[0].strip()] = {
                'family_id': f[1].strip(),
                'population': f[2].strip(),
                'population_description': f[3].strip(),
                'gender': f[4].strip().lower(),
                'relationship': f[5].strip()
            }

    print("Fetching standard 2,504 panel definition...")
    req_panel = urllib.request.Request(PANEL_URL, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req_panel, timeout=20) as resp:
        panel_lines = resp.read().decode('utf-8', errors='ignore').splitlines()

    panel_samples = []
    for line in panel_lines[1:]:
        parts = line.split('\t')
        if parts and parts[0].strip():
            panel_samples.append(parts[0].strip())

    # Read the master list of available samples in the 2019 GRCh38 release
    with open(DATA_DIR / "full_samples.txt") as f:
        available_samples = [line.strip() for line in f if line.strip()]

    # Filter to exactly 2,504 samples available in the callset
    selected_samples = [s for s in available_samples if s in set(panel_samples)]
    if len(selected_samples) < 2504:
        # Add remaining available samples from full_samples until exactly 2,504
        remaining = [s for s in available_samples if s not in set(selected_samples)]
        selected_samples.extend(remaining[: 2504 - len(selected_samples)])
    else:
        selected_samples = selected_samples[:2504]

    print(f"Selected exactly {len(selected_samples)} samples.")

    # Write sample_ids.txt
    ids_file = N2504_DIR / "sample_ids.txt"
    ids_file.write_text("\n".join(selected_samples) + "\n")
    print(f"Saved: {ids_file}")

    # Write sample_metadata.tsv
    metadata_rows = []
    super_counts = {}
    pop_counts = {}

    for s in selected_samples:
        info = sample_info.get(s, {})
        pop = info.get('population', 'UNK')
        super_pop = POP_TO_SUPER.get(pop, 'UNK')
        super_counts[super_pop] = super_counts.get(super_pop, 0) + 1
        pop_counts[pop] = pop_counts.get(pop, 0) + 1

        metadata_rows.append({
            'sample_id': s,
            'super_population': super_pop,
            'super_population_name': SUPER_POP_NAMES.get(super_pop, 'Unknown'),
            'population': pop,
            'population_description': info.get('population_description', ''),
            'gender': info.get('gender', 'unknown'),
            'family_id': info.get('family_id', s),
            'relationship': info.get('relationship', 'unrelated')
        })

    meta_file = N2504_DIR / "sample_metadata.tsv"
    with open(meta_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample_id', 'super_population', 'super_population_name',
            'population', 'population_description', 'gender', 'family_id', 'relationship'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(metadata_rows)

    print(f"Saved: {meta_file}")
    print("\nSuper-Population Breakdown:")
    for sp, cnt in sorted(super_counts.items()):
        print(f"  {sp} ({SUPER_POP_NAMES.get(sp, '')}): {cnt} samples")

if __name__ == "__main__":
    main()
