# Notes on changes to Jimmy's code

Task 3, extended from chromosomes 4 and 18 to all 22 autosomes. That takes the
off-target evaluation from 251 of the 2,887 validated CRISPOR candidates to
2,726. Two notebooks: seven genomes, and the full 2,504-individual panel.

## What changed in `utils/variant_utils.py`

- **Undetermined haplotypes are reported, not assumed to be reference.** The
  original returns reference sequence when a genotype is missing or unusable,
  which downstream counts as "no variant here". Each such case now raises a
  named flag and returns `None`.
- **Wider reconstruction window.** Variants are applied on a per-base grid over
  ±500 bp rather than spliced into the 23-mer, so deletions reaching in from
  outside the target are handled instead of silently dropped.
- **Two new off-target categories.** `SHIFTED` (mismatch count unchanged,
  position moved) and `UNCHANGED_NO_PAM` (sequence changed, no NGG either way)
  were both previously counted as `UNCHANGED`.
- **`INDETERMINATE` individual status**, so missing data is not reported as
  therapeutic escape; ancestry rates now divide by evaluable individuals.
- **Phasing check** that fails the notebook rather than producing
  haplotype-level results on unphased data.
- Reporting and plotting helpers for the genome-wide figures and tables.

Re-running both cohorts on this version reproduced every published result. The
silent-fallback path never fires on 1000 Genomes phased data, so the fixes are
insurance rather than corrections. The only visible differences are the
observations moved into the two new categories.

Unchanged from the `jimmy` branch: the VCF preparation scripts, the gene
architecture plotting, and the Task 2 inputs in `inputs/`.

## Running it

```bash
pip install -r requirements.txt
export IFN646_DATA=/path/to/data        # GRCh38.fa + vcfs/, too large for git
python3 utils/prepare_regional_vcfs.py --chromosomes autosomes
python3 utils/prepare_n2504_vcfs.py    --chromosomes autosomes
python3 utils/test_variant_utils.py     # 42 checks
```

`bcftools` is needed to build the VCFs, not to read the results.

## Scope

Autosomes only. The 1000 Genomes 20190312 release has no chromosome Y file, and
its chromosome X file covers only the pseudoautosomal regions, where none of the
candidates fall — so 161 of 2,887 sites cannot be evaluated. They are reported
as not evaluated, never as unchanged.
