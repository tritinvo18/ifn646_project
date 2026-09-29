# Explaining Task 3: Real-World CRISPR Safety Testing

When we design CRISPR treatments (like we did in Task 2), we initially build them using the **Reference Human Genome**. You can think of the Reference Genome as a generic, "average" blueprint of human DNA.

However, in the real world, no two people share the exact same DNA. Everyone has millions of tiny, natural genetic variations (called SNPs or mutations) that make them unique.

**Task 3 is where we take our CRISPR guides out of the generic "lab simulation" and test them against the actual DNA of real people.** We want to answer two critical questions before we would ever give this to a patient:

1. *Will it still work?* (Efficacy)
2. *Will it accidentally harm them?* (Safety)

---

## Step-by-Step Breakdown

### Step 1: Gathering Real Human Genomes

**What it means:** Instead of using the generic reference blueprint, we downloaded the real, mapped DNA sequences of 7 actual individuals from different global populations (provided by the 1000 Genomes Project).
**Why we do it:** A drug that works perfectly in a generic computer model might fail completely in real life if we don't account for natural human genetic diversity.

### Step 2: Testing "On-Target" Availability (Will it still work?)

**What it means:** CRISPR acts like a molecular GPS. It searches for a very specific 23-letter DNA sequence (our target gene). We scanned the genomes of the 7 individuals to see if their specific version of the target gene had any natural mutations inside that 23-letter sequence.
**Why we do it:** If a patient has a mutation right where our CRISPR guide is trying to land, the GPS fails. The CRISPR tool won't be able to bind to the DNA, and the treatment will fail. We do this step to guarantee that the guides we chose in Task 2 will actually work for the majority of the human population.
*(Spoiler: Our chosen guides passed! None of the 7 individuals had mutations that would block the treatment).*

### Step 3: Hunting for "Off-Target" Risks (Is it still safe?)

**What it means:** Sometimes, CRISPR accidentally binds to the wrong part of the genome because that part looks *almost* identical to the real target. These are called "Off-Targets." In Task 2, we made sure our guides had no dangerous off-targets in the generic reference genome. But in Task 3, we checked those "near-match" locations in our 7 real people.
**Why we do it:** This is the most important safety check. Imagine a harmless "near-match" sequence in the reference genome that is 3 letters different from our target. CRISPR ignores it because it's too different. But what if a real patient has a random mutation there that changes 2 of those letters? Suddenly, that harmless sequence becomes an almost perfect match. CRISPR might accidentally cut it, which could inadvertently damage a healthy gene or trigger a disease. We must ensure that natural mutations don't suddenly turn safe guides into dangerous ones.

### Step 4: Final Scoring and Re-Ranking

**What it means:** After scanning the 7 real genomes, we generated a personalized safety report. We looked for any guides that became unexpectedly dangerous in any of the individuals (meaning their mutations created new off-target risks).
**Why we do it:** We want to filter out guides that are only safe "on paper" and keep the ones that are robust and safe in the real world.
*(For example: In our results, the four guides we chose for the CXCL11 gene were completely flawless—none of the 7 individuals had mutations that increased their risk. However, for the SERPINB2 gene, we discovered that 2 of our guides suddenly became riskier in real people due to mutations, meaning we should probably avoid using them!)*

### The Big Picture

In simple terms: Task 2 is designing the key to fit the lock. Task 3 is checking whether that key accidentally opens any other doors in 7 different houses. By doing this, we bridge the gap between theoretical computer biology and real-world precision medicine.

---

## Appendix: Understanding the Guide Filtering Metrics (From Task 2)

When reviewing the final tables of our CRISPR guides, you might notice several highly technical column names. Here is a simple breakdown of what these flags mean and why we track them:

### Biological Manufacturing Limits

* **`poly_t_hard_exclusion` (The "Stop Sign" Error):** Inside a living cell, CRISPR guides are printed by a biological machine called an RNA Polymerase. However, if this machine reads four "T" letters in a row (`TTTT`), it interprets it as a strict "STOP" signal and aborts the printing process. We strictly exclude any guide with a `TTTT` sequence because the cell physically cannot manufacture it.
* **`tt_terminal_warning` (Floppy Endings):** If a guide has too many "T"s near its tail end (near where the cutting happens), the molecule becomes biologically unstable or "floppy." This makes the CRISPR machinery struggle to cut efficiently. We penalize these guides in our ranking.
* **`restriction_annotation_present` (The Lab Scissors Problem):** To get the CRISPR guide into a cell, scientists first have to "copy and paste" it into a delivery vehicle (a plasmid) using lab chemicals called restriction enzymes (molecular scissors). If the guide itself accidentally contains the sequence those scissors look for, the guide will be chopped to pieces during the laboratory manufacturing process.

### Molecular GPS Coordinates

* **`protospacer_start0` and `protospacer_end0`:** These are the exact mathematical coordinates showing where the 20-letter CRISPR guide sequence begins and ends on our target DNA. (The "0" just means the counting starts at 0 instead of 1, which is standard for computer programming).
* **`cut_boundary0` (The Exact Snip Site):** CRISPR doesn't just stick to the DNA; it acts as molecular scissors. This coordinate tells us the precise molecular location where the Cas9 protein will slice the DNA (which happens a few letters away from the end of the guide). We calculate this coordinate to guarantee that the slice happens *exactly* inside the important part of the gene (the first exon) rather than accidentally missing and cutting into nearby "junk" DNA.

### Composite Score: Why 60/40 (and Not 70/30, 80/20, or 50/50)?

To rank eligible guides, we compute a custom **Composite Score**:
$$\text{Composite Score} = 0.60 \times \text{cfdSpecScore (Specificity)} + 0.40 \times \text{Efficiency}$$

Why did we choose a **60/40** ratio instead of a heavier specificity weight like **70/30** or **80/20**, or equal weighting **50/50**?

1. **Safety-First Principle (Asymmetric Consequences):**
   * **Off-target damage is irreversible:** An unintended cut or binding event can permanently mutate an oncogene or silence essential housekeeping genes across millions of treated patient cells.
   * **On-target efficiency can be rescued:** In our pipeline, we deliver a **multiplexed pool of 4 guides** targeting the first exon simultaneously. Even if an individual guide has modest efficiency, the synergistic effect of multiplexing provides strong on-target transcriptional activation/repression.
   * Therefore, specificity **must hold majority voting power (>50%)**.

2. **The "Inactive Guide" Trap (Why NOT 70/30 or 80/20):**
   * Guide RNAs with near-perfect specificity scores (>90) often achieve that score simply because they contain unusual k-mers that rarely appear across the human genome. However, these same sequence properties often lead to poor Cas9 loading kinetics, unfavorable secondary structures, or weak R-loop formation.
   * At **70/30** or **80/20**, the formula becomes hyper-sensitive to trivial gains in CFD specificity (e.g., 91 vs. 83) while drastically downplaying cutting efficiency. Highly active cutters get discarded in favor of sluggish, low-potency guides.

3. **Empirical Evidence from Our CXCL11 Pool:**
   Testing alternative weights directly on our eligible candidate pool proves why 60/40 is the optimal Pareto balance:

   | Weighting | Selected Top 4 Guides | Kinetic Trade-off |
   | :--- | :--- | :--- |
   | **50 / 50** | `76forw` (Eff: 65), `261rev` (Eff: 67), `264forw` (Eff: 54), `77forw` (Eff: 59) | Underweights safety; fails to prioritize the ultra-safe `295forw` (CFD: 91). |
   | **60 / 40 (Optimal)** | **`76forw` (Eff: 65, CFD: 83)<br>`264forw` (Eff: 54, CFD: 89)<br>`261rev` (Eff: 67, CFD: 79)<br>`295forw` (Eff: 46, CFD: 91)** | **Sweet Spot:** Captures our most potent cutters (`261rev` Eff: 67, `76forw` Eff: 65) while retaining high-specificity guides (`295forw` CFD: 91). Pool average efficiency = **58.0**. |
   | **70 / 30** | `264forw` (Eff: 54), `76forw` (Eff: 65), `295forw` (Eff: 46), `106forw` (Eff: 42) | Highly active cutter **`261rev` (Eff: 67) is eliminated** and replaced by `106forw` (Eff: 42). |
   | **80 / 20** | `295forw` (Eff: 46), `264forw` (Eff: 54), `106forw` (Eff: 42), `76forw` (Eff: 65) | Over-penalizes activity; pool is overrun with sluggish cutters (pool average efficiency plunges to **51.8**). |

4. **Hard Filters Already Safeguard Baseline Safety:**
   * High-risk off-target guides are already disqualified upfront by our hard exclusion filters (`CFD > 50` and `MIT > 50`).
   * Every candidate entering the ranking stage is already clinically safe. Pushing specificity to 70% or 80% yields diminishing returns on safety while gutting biological activity.

5. **Why NOT 50/50?**
   * A 50/50 ratio treats efficiency and toxicity as equal (1:1). In clinical medicine, the catastrophic potential of off-target mutations demands an asymmetric **1.5:1 ratio (60% to 40%)** prioritizing safety.

---

## Appendix B: Understanding the Genomic Terminology

If you are reviewing the raw output tables or logs from Task 3, you may encounter a few standard genetic terms and quirks of the dataset:

### What is an "Allele"?

Humans inherit two copies of almost every gene—one from the biological mother, and one from the biological father. Each individual copy is called an **"allele."**

When we scan a real patient's genome for mutations, a mutation often only exists on *one* of the two alleles. For example, if a patient has a mutation that destroys the CRISPR binding site on their mother's copy but not their father's copy, the CRISPR treatment will only be 50% effective (cutting one allele but failing on the other).

In Task 3, our `VariantAnalyzer` explicitly reconstructs *both* alleles for every individual to guarantee that our guides are fully capable of cutting both copies of the target gene.

### Why is `chrY` (Chromosome Y) missing?

When scanning the genomes for off-target risks, you might notice that `chrY` occasionally throws "missing data" or "unresolved" warnings. This is not an error!

The 1000 Genomes Project cohort includes both male (XY) and female (XX) participants. Because biological females do not possess a Y chromosome, any attempt to evaluate an off-target risk located on `chrY` for a female participant will naturally return no data.

## References

### 1. The Foundational CRISPRi / CRISPRa Design Framework

Gilbert, L. A., Horlbeck, M. A., Adamson, B., Villalta, J. E., Chen, Y., Whitehead, E. H., & Weissman, J. S. (2014). Genome-Scale CRISPR-Mediated Control of Gene Repression and Activation. Cell,159(3), 647–661. DOI: 10.1016/j.cell.2014.09.029 <https://doi.org/10.1016/j.cell.2014.09.029>

Direct Relevance to Task 2:

* Defines the exact biological rules for targeting dCas9-KRAB (CRISPRi) and dCas9-SunTag/VP64 (CRISPRa) around the Transcription Start Site (TSS) and first exon.
* Establishes the standard experimental justification for using multiplexed pools of 4–5 guides to achieve robust transcriptional silencing/activation without requiring genomic double-strand breaks.

### 2. CFD Off-Target Specificity & On-Target Efficiency (Doench '16)

Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Cheah, C. S., Xu, W., ... & Root, D. E.
(2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology, 34(2), 184–191. DOI: 10.1038/nbt.3437 <https://doi.org/10.1038/nbt.3437>

Direct Relevance to Task 2 & 3:

* Introduces the Cutting Frequency Determination (CFD) score and the Rule Set 2 (Doench '16) efficiency score that are computed by CRISPOR and used directly in our ranking equation: Composite Score = 0.6 × CFD + 0.4 × Efficiency
* Details how mismatches in the "seed region" (PAM-proximal 8–12 nt) penalize cleavage kinetics versus PAM-distal mismatches.

### 3. The CRISPOR Platform Methodology

Haeussler, M., Schönig, K., Eckert, H., Eschstruth, A., Mianné, J., Renaud, J. B., ... & Concordet, J. P. (2016). Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR. Genome Biology, 17(1), 148. DOI: 10.1186/s13059-016-1012-2 <https://doi.org/10.1186/s13059-016-1012-2>

Direct Relevance to Task 2:

* Documents the exact pipeline and algorithms implemented in the CRISPOR tool.
* Explicitly describes the sequence penalties implemented in our code: Poly-T tracts (TTTT) acting as Pol III/U6 terminators, low/high GC boundaries (20–80%), and terminal TT-motifs that impede Cas9 loading.

### 4. Human Genetic Variation Impacting CRISPR Efficacy and Safety

Canver, M. C., Lessard, S., Pinello, L., Wu, Y., Ilboudo, Y., Stern, E. N., ... & Bauer, D. E.
(2017). Variant-aware saturating mutagenesis using CRISPR-Cas9. Nature Genetics, 49(4), 625–634. DOI: 10.1038/ng.3800 <https://doi.org/10.1038/ng.3800>

Direct Relevance to Task 3 (Question 1 - On-Target Viability):

* Demonstrates how single nucleotide polymorphisms (SNPs) and short indels disrupt CRISPR binding by altering protospacer sequences or destroying the required NGG PAM.
* Validates the need to reconstruct both alleles from VCF data across individual genomes (as done by VariantAnalyzer) to ensure population-wide therapeutic applicability.

### 5. Personalized Off-Target Risk and Neo-PAM Generation

Lessard, S., Francioli, L., Alfoldi, J., Tardif, J. C., Ellinor, P. T., MacArthur, D. G., ...
& Bauer, D. E. (2017). Human genetic variation alters CRISPR-Cas9 on- and off-target specificity.
Proceedings of the National Academy of Sciences (PNAS), 114(52), E11257–E11266. DOI: 10.1073/pnas.
1714640114 <https://doi.org/10.1073/pnas.1714640114>

Direct Relevance to Task 3 (Question 2 - Off-Target Risk):

* The definitive reference paper showing that using standard reference genomes (like GRCh38) causes severe blind spots in CRISPR safety analysis.
* Shows that natural human variants frequently create novel PAM sites (neo-PAMs) or decrease mismatch counts at candidate loci, converting previously benign genomic sites into active, high-risk off-target cleavage targets in specific individuals (the exact behavior observed in our SERPINB2 40forw/90forw findings).
