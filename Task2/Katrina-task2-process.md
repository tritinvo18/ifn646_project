# SERPINB2 exon 1 — sequence extraction


---

## Result

| | |
|---|---|
| Gene | SERPINB2 |
| Assembly | GRCh38 (`GCF_000001405.26`) |
| Accession | `NC_000018.10` (chromosome 18) |
| Strand | `+` |
| Gene extent | 63,887,705 – 63,903,890 |
| Exon 1 | 63,887,705 – 63,887,770 (66 bp) |
| Extracted window | 63,887,680 – 63,887,795 (exon ± 25 bp) |

```
AACCAGTCATTACCATGTCTGAACTGTAACAACTCTCAGAGGAGCATTGCCCGTCAGACAGCAACTCAGAGAATAACCAGAGAACAACCAGGTATTTCAATGATTTCCATGCCATG
```

Within this 116 bp window: positions 1–25 are upstream, 26–91 are exon 1,
92–116 are intron 1. Position 26 is the TSS. 

Note - I think this differs from Jonathan's window. 

With this result I used the CHOPCHOP web interface. Results are saved in a separate folder in my branch. 

---

## Findings

**Both annotated transcripts share an identical first exon.** `NM_001143818.1`
and `NM_002575.2` span the same extent and both begin at 63,887,705 –
63,887,770. They differ only by an additional internal exon
(63,889,838 – 63,890,095) present in variant 1. The transcript choice therefore
has no bearing on this task.

No MANE Select tag is available: `GCF_000001405.26` is the original GRCh38
release, which predates MANE. The transcripts were compared directly instead.

**25 bp of flanking sequence was included either side.** A 20 nt protospacer
overlapping the exon may begin up to 19 bases before it, so the search space is
the exon plus flank rather than the 66 bp alone. This matters here because the
exon is short.

I checked the strand, coordinates and exon 1 boundaries were cross-checked against NCBI Gene via the website.

---

## Bash commands and process

Genome files from the NCBI FTP site listed in the project brief
(`GCF_000001405.26_GRCh38_genomic.fna.gz` and `.gff.gz`). Sequence extraction
used **pyfaidx** rather than samtools; equivalent functionality, index built
automatically on first use.

Locate the gene and read its strand:

```bash
gunzip -c GCF_000001405.26_GRCh38_genomic.gff.gz \
  | awk -F'\t' '$3=="gene" && $9 ~ /Name=SERPINB2;/'
```

List transcripts:

```bash
gunzip -c GCF_000001405.26_GRCh38_genomic.gff.gz \
  | awk -F'\t' '$3=="mRNA" && $9 ~ /gene=SERPINB2/' \
  | cut -f4,5,9
```

Compare exon structures across both transcripts:

```bash
for tx in NM_001143818.1 NM_002575.2; do
  echo "== $tx"
  gunzip -c GCF_000001405.26_GRCh38_genomic.gff.gz \
    | awk -F'\t' -v t="$tx" '$3=="exon" && $9 ~ ("Parent=rna-" t ";")' \
    | cut -f4,5 | sort -k1,1n
done
```

Confirm where translation starts:

```bash
gunzip -c GCF_000001405.26_GRCh38_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS" && $9 ~ /Parent=rna-NM_002575.2;/' \
  | cut -f4,5 | sort -k1,1n | head -3
```

Extract the sequence:

```bash
gunzip GCF_000001405.26_GRCh38_genomic.fna.gz
faidx GCF_000001405.26_GRCh38_genomic.fna NC_000018.10:63887680-63887795 \
  > serpinb2_exon1_ctx.fa
```

Verify (expect 116 bases, no ambiguous positions):

```bash
grep -v ">" serpinb2_exon1_ctx.fa | tr -d '\n' | wc -c
grep -v ">" serpinb2_exon1_ctx.fa | grep -c "N"
```

GFF and faidx both use 1-based inclusive coordinates, so GFF values were used
unchanged. The gene is on the plus strand, so no reverse-complementing was
required.

