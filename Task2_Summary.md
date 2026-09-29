# Task 2 Summary

## What is a "Hard Filter" (eligible = False)?

A guide is only explicitly disqualified (marked eligible = False) if it violates an absolute
physical or biological constraint:

- Location: The Cas9 cut site is physically outside the first exon boundary.
- Safety: It fails the minimum safety threshold (CFD/MIT < 50), risking lethal off-target cleavage.
- Transcription Failure: It contains a Poly-T tract (TTTT). Note: This is the only "Inefficient" flag that acts as a hard filter, because U6/U3 promoters will physically stop transcribing the RNA when they hit TTTT, resulting in a broken, non-functional guide.

### How are "Inefficient" guides handled? (Soft Filters)

Guides that CRISPOR flags as merely "Inefficient" (due to bad GC content, TT-rich regions near the PAM, or low Doench '16 scores) are kept in the eligible pool, but mathematically penalized in the ranking algorithm:

They are ranked using this logic:

```python
ranked = guides.sort_values(
    ["eligible", "tt_terminal_warning", "Composite_Score", "cfdSpecScore", "Efficiency"],
    ascending=[False, True, False, False, False],
)
```

Because of this sorting order:

1. Guides with a TT-terminal warning are pushed down below any guides without one.
2. The Composite Score factors in the Doench '16 Efficiency score. Guides with poor efficiency
naturally sink to the bottom of the list.
