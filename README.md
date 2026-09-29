# IFN646 Project

This repository contains the analysis codebase for the IFN646 project, focusing on identifying dysregulated genes in COVID-19 and designing targeted CRISPR strategies.

## Task 3: Individual-Genome CRISPR Guide Evaluation

Task 3 requires evaluating our selected Task 2 CRISPR guides (for SERPINB2 and CXCL11) across multiple individual human genomes to ensure they remain usable and safe despite naturally occurring genetic variation.

### Downloading Required Data

The project `REQUIREMENTs.md` states that preprocessed variants for seven individuals are provided for Chromosome 1 and Chromosome 2 (`sampled_chr1.vcf` and `sampled_chr2.vcf`).

However, our selected target genes reside on different chromosomes:

* **CXCL11** is located on **Chromosome 4**
* **SERPINB2** is located on **Chromosome 18**

Additionally, the CRISPOR off-target candidates predicted in Task 2 are distributed across all autosomes and the X chromosome. Therefore, to faithfully reproduce the Task 3 off-target and on-target validation, we must obtain variant data for these specific seven individuals across all relevant chromosomes.

#### 1. Extracting Target Sample IDs

First, extract the seven individual sample IDs directly from the provided subset files:

```bash
bcftools query -l data/sampled_chr1.vcf > data/task3_sample_ids.txt
```

#### 2. Automated Regional VCF Downloading

Instead of manually downloading the entire 1000 Genomes VCFs for every chromosome (which are enormous), we have provided a utility script that streams and downloads only the necessary target and off-target regions for our 7 specific individuals directly from the 1000 Genomes FTP server.

Run the preparation script from the project root:

```bash
python utils/prepare_regional_vcfs.py --chromosomes autosomes
```

This script will:

1. Connect to the remote 1000 Genomes FTP server.
2. Query the regions intersecting our guides' on-target locations (chr4, chr18) and all predicted off-target sites.
3. Filter the remote data to extract **only** our seven target samples.
4. Normalize the variants and enforce reference sequence conformity against our local GRCh38 FASTA (`data/GRCh38.fa`).
5. Save the final output as compact bgzipped and CSI-indexed VCF files in the `data/vcfs/` directory (e.g., `data/vcfs/sampled_chr4.norm.vcf.gz`).

### Reproducing the Task 3 Analysis

Once the required chromosome VCFs are populated in `data/vcfs/`, you can reproduce the complete Task 3 analysis:

1. **Setup the Environment:**
   Ensure you have the required Python packages installed:

   ```bash
   pip install pysam pandas numpy matplotlib seaborn xlrd notebook
   ```

2. **Execute the Notebook:**
   The entire analysis is self-contained within the `03_task3_scientific_report.ipynb` Jupyter Notebook.
   You can run it interactively by starting Jupyter:

   ```bash
   jupyter notebook notebooks/03_task3_scientific_report.ipynb
   ```

   Or execute it directly from the command line to regenerate the final markdown report:

   ```bash
   jupyter nbconvert --to markdown --no-input --execute notebooks/03_task3_scientific_report.ipynb --output task3_report.md
   ```

3. **Outputs:**
   * Detailed TSV tables representing guide validation, usability, and off-target risk changes will be deposited in `results/tables/`.
   * Resulting distribution heatmaps will be deposited in `results/figures/`.
   * The final interpreted text report and per-guide summary will be generated as `notebooks/task3_report.md`.
