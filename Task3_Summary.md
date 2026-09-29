# Summary

Based on the experiment setup in [Task 3](/notebooks/03_task3_scientific_report.ipynb),both CXCL11 and SERPINB2 were successfully taken forward for guide-level evaluation. All eight selected guides mapped uniquely to GRCh38, retained compatible PAMs, and remained targetable on both alleles in all seven individuals. Therefore, neither CXCL11 CRISPRi nor SERPINB2 CRISPRa was limited by loss of the intended target in this cohort. However, individual genetic variation altered the predicted sequence-level off-target profiles of several guides.

The CXCL11 guide pool showed a flawless personalised off-target profile across this cohort. Zero CXCL11 sample-guide combinations contained strengthened candidate sites; all four guides (76forw, 264forw, 261rev, 295forw) were completely unaffected by risk-increasing genetic variation. By comparison, SERPINB2 40forw and 90forw showed strengthened candidate sites across all seven individuals, indicating that these guides should be deprioritised or subjected to additional validation. SERPINB2 50rev and 77rev had the most favourable profiles within the SERPINB2 pool.

These results support CXCL11 as the stronger candidate for subsequent functional investigation, as its entire 4-guide pool is highly robust against off-target variance. SERPINB2 remains technically targetable, but its biological role is less directly supported and its original four-guide pool should be refined, with 50rev and 77rev retained as the leading guides and 40forw and 90forw treated cautiously. Nevertheless, these results describe predicted guide suitability rather than therapeutic effectiveness. Functional experiments are required to determine whether CXCL11 inhibition or SERPINB2 activation produces the intended biological effect without impairing antiviral responses or causing off-target gene regulation

## Both genes can be targeted in all seven individuals

All four CXCL11 CRISPRi guides and all four SERPINB2 CRISPRa guides:

- mapped uniquely to GRCh38,
- had compatible NGG PAMs,
- retained both targetable alleles,
- and remained usable in all seven analysed genomes.

Therefore, none of the seven individuals would be excluded from either strategy because of an on-target sequence variant. Both four-guide pools retained complete predicted on-target coverage.

However, this result applies only to these seven genomes. It does not establish universal population coverage.

## Genetic variation affected predicted off-target risk, but not on-target availability

The analysis evaluated 1,457 of the 1,541 CRISPOR candidate off-target sites, corresponding to 94.55% of the selected candidate set. Across these sites, some individual variants strengthened potential off-target matches, while others weakened or removed them. Overall, created-or-strengthened allele observations occurred across 14 sample-guide combinations (all within SERPINB2). No completely new off-target sites were created among the evaluated candidates, but some existing candidate sites became more similar to their guide sequences.

Therefore, individual genetic variation did not alter whether the guides could bind their intended targets, but it did produce individual-specific differences in predicted off-target risk.

This directly answers the two genomic questions:

1. Were any guides absent or unusable? No.
2. Did predicted off-target risk differ among genomes? Yes.

## Comparison between the two gene strategies

### CXCL11 CRISPRi showed an incredibly consistent off-target profile

Remarkably, **zero** CXCL11 sample-guide combinations contained a strengthened candidate site. All four guides in the finalized CXCL11 pool (76forw, 264forw, 261rev, 295forw) returned `NO_SEQUENCE_LEVEL_CHANGE_DETECTED` across all seven evaluated genomes.

In contrast to early iterations of the pool, the strict enforcement of first-exon boundaries during Task 2 yielded an exceptionally robust set of guides. None of the four guides suffered from individual genetic variation increasing similarity at known off-target loci.

Thus, within the CXCL11 pool:

All four guides (76forw, 264forw, 261rev, 295forw) appear to have perfectly clean personalised sequence-level profiles across this cohort.

The important caveat is that the reported maximum CFD values apply to the reference candidates. They were not recalculated as personalised CFD scores. Consequently, guide ranking should consider both reference CRISPOR results and the personalised sequence changes, rather than treating one metric as definitive.

### SERPINB2 CRISPRa showed more frequent personalised changes

All of the sample-guide combinations with strengthened candidates in this experiment belonged to the SERPINB2 pool. In particular:

SERPINB2 40forw had strengthened candidates in all seven individuals.
SERPINB2 90forw had strengthened candidates in all seven individuals.
SERPINB2 50rev had no strengthened candidates in any individual.
SERPINB2 77rev had no strengthened candidates in any individual.

This indicates that the elevated SERPINB2 pool-level signal is entirely driven by 40forw and 90forw, rather than all four SERPINB2 guides.

Within the SERPINB2 pool:

50rev and 77rev appear to have the cleanest personalised sequence-level profiles.
40forw and 90forw should be deprioritised or individually reviewed, because strengthened candidate sites (increased off-target risk) were detected consistently across the cohort.

This does not prove that 40forw or 90forw will generate harmful off-target activity. It means they have less favourable computational profiles and warrant more scrutiny before experimental use.

## Does this mean CXCL11 should be chosen over SERPINB2?

From a guide-design perspective: CXCL11 is currently favoured because it has:

- complete on-target retention,
- fewer strengthened off-target candidates,
- fewer sample-guide combinations showing increased risk,
- and at least two guides with relatively stable personalised profiles.

The SERPINB2 set also retains complete on-target coverage, but 40forw and 90forw repeatedly show strengthened candidates across individuals. Therefore, the CXCL11 pool appears computationally more robust than the complete SERPINB2 pool.

However, from a biological perspective, the experiments do not decide between them because they assessed:

- guide mapping,
- target-sequence conservation,
- PAM retention,
- candidate off-target sequence similarity,
- and variation among seven individual genomes.

Whereas:

- whether CXCL11 inhibition reduces harmful inflammation,
- whether it interferes with antiviral immunity,
- whether SERPINB2 activation is protective or harmful,
- whether either intervention changes viral replication,
- whether either intervention improves a cellular phenotype,
- or whether the predicted off-target sites are actually bound or regulated.

were not assessed. Therefore, this experiment can only answer: CXCL11 currently has stronger disease-specific biological justification and a cleaner personalised guide profile, whereas SERPINB2 remains a more exploratory biological candidate with a less favourable four-guide pool.

## Recommendation

1. Primary experimental candidate: CXCL11 CRISPRi.
2. Leading CXCL11 guides: 76forw, 264forw, 261rev, 295forw (all highly robust).
3. Secondary exploratory candidate: SERPINB2 CRISPRa.
4. Leading SERPINB2 guides: 50rev and 77rev.
5. Guides requiring caution: SERPINB2 40forw and 90forw.

## Future direction

Measure target expression, functional phenotype, and expression at the highest-risk off-target genes.
