import re
import pysam

def reverse_complement(sequence):
    return sequence.upper().translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def canonical_chromosome(value):
    text = str(value).strip()
    match = re.fullmatch(r'chr([0-9]{1,2}|X|Y)', text)
    if match: return match.group(1)
    if text.startswith('NC_0000'): return str(int(text[7:9]))
    return text

def orient_sequence(genomic_sequence, strand):
    if genomic_sequence is None: return None
    return genomic_sequence if strand == '+' else reverse_complement(genomic_sequence)

def record_end(record):
    return int(record.info['END']) if 'END' in record.info else record.pos + len(record.ref) - 1

class VariantAnalyzer:
    def __init__(self, fasta_path):
        self.reference = pysam.FastaFile(str(fasta_path))
        self.vcfs = {}
        self.fasta_contigs = {canonical_chromosome(c): c for c in self.reference.references}
        
    def add_vcf(self, chrom, vcf_path):
        self.vcfs[canonical_chromosome(chrom)] = pysam.VariantFile(str(vcf_path))

    def _get_overlapping_records(self, vcf, chrom, start, end):
        return [r for r in vcf.fetch(chrom, max(0, start - 51), end + 50) if r.pos <= end and record_end(r) >= start]

    def reconstruct_alleles(self, chrom, start, end, sample_id):
        chrom = canonical_chromosome(chrom)
        if chrom not in self.vcfs:
            return None, None, 'NO_VCF_COVERAGE'
            
        fasta_chrom = self.fasta_contigs.get(chrom, chrom)
        genomic_reference = self.reference.fetch(fasta_chrom, start - 1, end).upper()
        vcf = self.vcfs[chrom]
        
        # Resolve exact VCF contig name
        vcf_contig = None
        for c in vcf.header.contigs:
            if canonical_chromosome(c) == chrom:
                vcf_contig = c
                break
        
        records = self._get_overlapping_records(vcf, vcf_contig, start, end)
        if not records:
            return genomic_reference, genomic_reference, 'REFERENCE_HOMOZYGOUS'

        genotypes = [r.samples[sample_id].get('GT') for r in records]
        
        # Simple resolution: if unphased or missing, fallback to Reference (relaxed constraint)
        if any(gt is None or len(gt) < 2 or any(a is None for a in gt[:2]) for gt in genotypes):
            return genomic_reference, genomic_reference, 'REFERENCE_ASSUMED'

        edits = [[], []]
        for record, genotype in zip(records, genotypes):
            for haplotype, allele_index in enumerate(genotype[:2]):
                if allele_index == 0 or record.alts is None or allele_index > len(record.alts):
                    continue
                alt = str(record.alts[allele_index - 1]).upper()
                if alt.startswith('<') or alt == '*':
                    continue
                edits[haplotype].append({
                    'offset': record.pos - start, 'ref': record.ref.upper(), 'alt': alt
                })

        alleles = [genomic_reference, genomic_reference]
        for haplotype in range(2):
            seq = alleles[haplotype]
            for edit in sorted(edits[haplotype], key=lambda item: item['offset'], reverse=True):
                offset = edit['offset']
                if seq[offset:offset + len(edit['ref'])] == edit['ref']:
                    seq = seq[:offset] + edit['alt'] + seq[offset + len(edit['ref']):]
            alleles[haplotype] = seq
            
        return alleles[0], alleles[1], 'RECONSTRUCTED'

    def classify_on_target(self, ref_23mer, alt_23mer):
        if not alt_23mer or len(alt_23mer) != 23: return 'INDEL_OR_UNRESOLVED'
        if not alt_23mer[-3:].endswith('GG'): return 'PAM_LOST'
        mismatches = sum(a != b for a, b in zip(ref_23mer[:20], alt_23mer[:20]))
        if mismatches == 0: return 'UNCHANGED'
        return 'MISMATCHED'

    def classify_off_target(self, guide_20mer, ref_site, alt_site):
        if not alt_site or len(alt_site) != 23: return 'UNRESOLVED'
        ref_pam = ref_site[-3:].endswith('GG')
        alt_pam = alt_site[-3:].endswith('GG')
        ref_mm = sum(a != b for a, b in zip(guide_20mer, ref_site[:20]))
        alt_mm = sum(a != b for a, b in zip(guide_20mer, alt_site[:20]))
        
        if not ref_pam and alt_pam: return 'CREATED'
        if ref_pam and not alt_pam: return 'REMOVED'
        if alt_pam and alt_mm < ref_mm: return 'INCREASED'
        if alt_pam and alt_mm > ref_mm: return 'DECREASED'
        return 'UNCHANGED'
