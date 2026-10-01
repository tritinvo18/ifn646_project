# Population-Scale CRISPR Variant Analysis Report (1000 Genomes N = 2,504, All 14 Eligible Guides)

## 1. Executive Summary

This study executes a population-scale genomic evaluation of candidate CRISPR interference (CRISPRi) and CRISPR activation (CRISPRa) guide RNAs across the complete **1000 Genomes Project Phase 3 cohort ($N = 2,504$ individuals, 5,008 haplotypes)**. Rather than relying on the pedagogical 7-sample toy panel or truncating candidates to an arbitrary pool size ($k=4$), this expanded analysis evaluates all **14 eligible gRNAs** identified in Task 2 (9 *CXCL11* CRISPRi guides and 5 *SERPINB2* CRISPRa guides) across five major continental super-populations: **African (AFR)**, **European (EUR)**, **East Asian (EAS)**, **South Asian (SAS)**, and **Admixed American (AMR)**.

### Primary Breakthroughs & Biological Insights:
1. **The East Asian Blindspot in *SERPINB2* (`50rev`)**:
   - In the initial 7-sample panel (which lacked East Asian representation), guide `50rev` appeared 100% functional.
   - In the $N = 2,504$ cohort, `50rev` exhibits **300 mismatched alleles**, heavily concentrated in East Asian donors (**33.93% of East Asian individuals disrupted**; allele frequency $> 20\%$). In an East Asian therapeutic trial, over one in three patients would suffer reduced or abolished CRISPRa activation.
2. **African-Specific PAM Loss in *SERPINB2* (`51rev`)**:
   - Guide `51rev` exhibits **16 canonical `NGG` PAM loss events**, 15 of which occur exclusively in African individuals (`AFR`), rendering them completely refractory to Cas9 binding.
3. **Identification of Universal, Ancestry-Robust Therapeutic Guides**:
   - On **Chromosome 4 (*CXCL11*, CRISPRi)**: Guides **`106forw`**, **`145rev`**, **`169rev`**, and **`295forw`** achieve **100.0% biallelic targetability** (2,504/2,504 individuals, 5,008/5,008 intact alleles across all 5 continental groups).
   - On **Chromosome 18 (*SERPINB2*, CRISPRa)**: Guides **`40forw`** and **`77rev`** achieve **100.0% biallelic targetability** across all 2,504 individuals with zero increased off-target risk.

---

## 2. Cohort Demographics and Ancestral Architecture

All 2,504 individuals were extracted from the 2019 GRCh38 re-aligned release using remote indexed `bcftools view -R` queries via [`utils/prepare_n2504_vcfs.py`](file:///home/jimmy/Project/ifn646_project/utils/prepare_n2504_vcfs.py). The demographic structure encompasses:

| Continental Super-Population | Population Code | Included Sub-Populations (26 Total) | Sample Count ($N$) | Haplotypes ($2N$) |
| :--- | :--- | :--- | :--- | :--- |
| **African** | `AFR` | YRI, LWK, GWD, MSL, ESN, ASW, ACB | **660** | 1,320 |
| **European** | `EUR` | CEU, TSI, FIN, GBR, IBS | **504** | 1,008 |
| **East Asian** | `EAS` | CHB, JPT, CHS, CDX, KHV | **504** | 1,008 |
| **South Asian** | `SAS` | GIH, PJL, BEB, STU, ITU | **489** | 978 |
| **Admixed American** | `AMR` | MXL, PUR, CLM, PEL | **347** | 694 |
| **TOTAL** | — | **26 Global Populations** | **2,504** | **5,008** |

In contrast to the baseline 7-sample set (which contained 4 Europeans, 3 Africans, and 0 East Asians, South Asians, or Admixed Americans), the $N = 2,504$ cohort possesses the statistical power to detect low-frequency variants down to $MAF \approx 0.02\%$.

---

## 3. On-Target Variant Sensitivity & Biallelic Penetrance

Across $14 \text{ guides} \times 2,504 \text{ samples} \times 2 \text{ haplotypes} = \mathbf{70,112 \text{ reconstructed alleles}}$, we classified every target allele as `UNCHANGED`, `MISMATCHED`, `PAM_LOST`, or `INDEL_OR_UNRESOLVED`. Each individual's therapeutic status was scored as:
- **`PASS` (Biallelic Intact)**: 2 targetable alleles (100% theoretical efficacy).
- **`CAUTION` (Heterozygous Disrupted)**: 1 targetable allele (partial ~50% knockdown/activation).
- **`FAIL` (Homozygous Escape)**: 0 targetable alleles (complete therapeutic resistance).

### Comprehensive Biallelic Penetrance Ranking:

| Gene | Mode | Guide ID | Cohort Size ($N$) | Biallelic Pass ($N$) | Biallelic Pass (%) | Heterozygous Caution (%) | Homozygous Fail (%) | PAM Lost Alleles | Mismatched Alleles | Indel Alleles | Universal Penetrance Rank |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CXCL11** | CRISPRi | **`106forw`** | 2,504 | 2,504 | **100.0%** | 0.0% | 0.0% | 0 | 0 | 0 | **1 (Optimal)** |
| **CXCL11** | CRISPRi | **`145rev`** | 2,504 | 2,504 | **100.0%** | 0.0% | 0.0% | 0 | 0 | 0 | **1 (Optimal)** |
| **CXCL11** | CRISPRi | **`169rev`** | 2,504 | 2,504 | **100.0%** | 0.0% | 0.0% | 0 | 0 | 0 | **1 (Optimal)** |
| **CXCL11** | CRISPRi | **`295forw`** | 2,504 | 2,504 | **100.0%** | 0.0% | 0.0% | 0 | 0 | 0 | **1 (Optimal)** |
| **CXCL11** | CRISPRi | `261rev` | 2,504 | 2,503 | 99.96% | 0.04% | 0.0% | 1 | 0 | 0 | 5 |
| **CXCL11** | CRISPRi | `264forw` | 2,504 | 2,503 | 99.96% | 0.04% | 0.0% | 0 | 1 | 0 | 5 |
| **CXCL11** | CRISPRi | `98forw` | 2,504 | 2,502 | 99.92% | 0.08% | 0.0% | 0 | 2 | 0 | 7 |
| **CXCL11** | CRISPRi | `76forw` | 2,504 | 2,501 | 99.88% | 0.12% | 0.0% | 0 | 2 | 1 | 8 |
| **CXCL11** | CRISPRi | `77forw` | 2,504 | 2,501 | 99.88% | 0.12% | 0.0% | 0 | 2 | 1 | 8 |
| **SERPINB2**| CRISPRa | **`40forw`** | 2,504 | 2,504 | **100.0%** | 0.0% | 0.0% | 0 | 0 | 0 | **1 (Optimal)** |
| **SERPINB2**| CRISPRa | **`77rev`** | 2,504 | 2,504 | **100.0%** | 0.0% | 0.0% | 0 | 0 | 0 | **1 (Optimal)** |
| **SERPINB2**| CRISPRa | `90forw` | 2,504 | 2,494 | 99.60% | 0.40% | 0.0% | 0 | 10 | 0 | 3 |
| **SERPINB2**| CRISPRa | `51rev` | 2,504 | 2,488 | 99.36% | 0.64% | 0.0% | 16 | 0 | 0 | 4 |
| **SERPINB2**| CRISPRa | `50rev` | 2,504 | 2,228 | **88.98%** | **10.06%** | **0.96%** | 0 | 300 | 0 | **5 (High Risk)** |

---

## 4. Ancestry Stratification and Population-Specific Failure Modes

Disruption rates were computed separately across the five continental super-populations ([`results/tables/n2504/task3_n2504_ancestry_stratification.tsv`](file:///home/jimmy/Project/ifn646_project/results/tables/n2504/task3_n2504_ancestry_stratification.tsv)):

```
                        Ancestry Disruption Rate (%)
Guide ID   Gene       European (EUR)  African (AFR)  East Asian (EAS)  South Asian (SAS)  Admixed Amer. (AMR)
106forw    CXCL11         0.00%           0.00%           0.00%             0.00%               0.00%
145rev     CXCL11         0.00%           0.00%           0.00%             0.00%               0.00%
169rev     CXCL11         0.00%           0.00%           0.00%             0.00%               0.00%
295forw    CXCL11         0.00%           0.00%           0.00%             0.00%               0.00%
261rev     CXCL11         0.00%           0.00%           0.20%             0.00%               0.00%
264forw    CXCL11         0.00%           0.00%           0.20%             0.00%               0.00%
76forw     CXCL11         0.00%           0.00%           0.40%             0.00%               0.29%
77forw     CXCL11         0.00%           0.00%           0.40%             0.00%               0.29%
98forw     CXCL11         0.00%           0.30%           0.00%             0.00%               0.00%
40forw     SERPINB2       0.00%           0.00%           0.00%             0.00%               0.00%
77rev      SERPINB2       0.00%           0.00%           0.00%             0.00%               0.00%
90forw     SERPINB2       0.00%           0.00%           1.98%             0.00%               0.00%
51rev      SERPINB2       0.00%           2.27%           0.00%             0.00%               0.29%
50rev      SERPINB2       4.76%           0.45%          33.93%             9.82%               8.65%
```

### Critical Ancestry Observations:
1. **The East Asian Mismatch Hotspot (`50rev`)**:
   - `50rev` carries a high-frequency single-nucleotide variant in the protospacer region that is common in East Asia ($MAF \approx 20\%$) but rare in Africa ($MAF < 0.5\%$).
   - In East Asians, **33.93% of patients (171/504)** carry at least one disrupted allele, and **24 patients** are homozygous escapees (`FAIL`).
   - Screening only European and African individuals masked this defect entirely.
2. **The African PAM Loss Hotspot (`51rev`)**:
   - `51rev` carries an African-specific single-nucleotide variant that mutates the canonical `NGG` PAM to `NGA`/`NGT`.
   - 15 out of 16 observed PAM loss events occur in African donors (2.27% disruption rate in AFR; 0.00% in EUR, EAS, SAS).

---

## 5. Candidate Off-Target Cleavage Alterations

Across 251 candidate off-target loci on Chromosomes 4 and 18, evaluating all 2,504 individuals yielded **1,257,008 pairwise sequence comparisons**:
- **Unchanged Reference Sites**: $98.36\%$ (1,236,395 evaluations).
- **Altered Cleavage Events**: $1.64\%$ (20,613 events saved to [`results/tables/n2504/task3_n2504_personalised_offtargets.tsv`](file:///home/jimmy/Project/ifn646_project/results/tables/n2504/task3_n2504_personalised_offtargets.tsv)):
  - **`DECREASED` Risk**: 11,311 events (variants increase mismatch distance or reduce homology).
  - **`INCREASED` Risk**: 8,404 events (variants eliminate a mismatch, increasing similarity to the gRNA).
  - **`REMOVED` Sites**: 88 events (variants abolish the `NGG` PAM at the off-target locus, eliminating the cleavage hazard).
  - **`UNRESOLVED` Sites**: 810 events.

While the 7-sample analysis detected zero altered off-targets, the $N = 2,504$ analysis confirms that human variants frequently modify off-target affinities, reinforcing the need to prioritize guides with zero or minimal increased-risk events.

---

## 6. Visual Exhibit Gallery

All figures were generated at 300 DPI and nested cleanly in [`results/figures/n2504/`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/):

1. **Gene Sequence Architecture**:
   - *CXCL11* (Chr 4, 9 guides, Exon 1, PAM): [`task3_n2504_cxcl11_grna_architecture.png`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/task3_n2504_cxcl11_grna_architecture.png)
   - *SERPINB2* (Chr 18, 5 guides, Exon 1, PAM): [`task3_n2504_serpinb2_grna_architecture.png`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/task3_n2504_serpinb2_grna_architecture.png)
2. **Ancestry Disruption Heatmap**:
   - Matrix showing gRNA disruption rates across EUR, AFR, EAS, SAS, and AMR: [`task3_n2504_ancestry_disruption_heatmap.png`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/task3_n2504_ancestry_disruption_heatmap.png)
3. **Biallelic Targetability Waterfall**:
   - Stacked horizontal bar chart showing Biallelic Intact, Heterozygous Disrupted, and Homozygous Escape proportions: [`task3_n2504_biallelic_penetrance_waterfall.png`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/task3_n2504_biallelic_penetrance_waterfall.png)
4. **Apparent vs. True Targetability ($N=7$ vs. $N=2,504$)**:
   - Juxtaposition illustrating how small cohorts mask major population defects: [`task3_n2504_cohort_comparison_n7_vs_n2504.png`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/task3_n2504_cohort_comparison_n7_vs_n2504.png)
5. **Off-Target Effect Distribution**:
   - Classification breakdown of the 20,613 altered off-target events: [`task3_n2504_offtarget_risk_breakdown.png`](file:///home/jimmy/Project/ifn646_project/results/figures/n2504/task3_n2504_offtarget_risk_breakdown.png)

---

## 7. Definitive Clinical Recommendations

When designing CRISPR therapeutics for human viral hyperinflammation (such as COVID-19 in GSE159717), guide selection must guarantee **pan-ethnic efficacy** and **zero genotoxicity**:

### Recommendation 1: *CXCL11* Knockdown (CRISPRi)
- **Top Choice**: **`106forw`** (Rank 1 in Task 2 Composite Score: 74.4; 100.0% universal biallelic targetability across all 2,504 individuals; zero PAM disruptions; zero increased-risk off-targets).
- **Alternative Universal Options**: **`145rev`**, **`169rev`**, and **`295forw`** are equally 100.0% intact across all 5 continental groups and provide excellent multiplexing backups.
- **De-prioritized Guides**: `76forw` and `77forw` (carry indels in Admixed Americans and mismatches in East Asians); `261rev` (carries a PAM loss variant).

### Recommendation 2: *SERPINB2* Activation (CRISPRa)
- **Top Choice**: **`40forw`** (Top Task 2 Composite Score: 78.4; 100.0% universal biallelic targetability across all 2,504 individuals; zero on-target disruptions in any ancestry; safe off-target profile).
- **Alternative Universal Option**: **`77rev`** (100.0% biallelic targetability across all 2,504 individuals).
- **STRICTLY EXCLUDED GUIDES**:
  - **`50rev`**: Rejected due to severe East Asian failure mode (33.93% disruption rate, 24 homozygous escapees).
  - **`51rev`**: Rejected due to African-specific PAM loss (15 PAM mutations in AFR donors).
