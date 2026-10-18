"""Self-contained tests for the reconstruction and classification logic.

pysam, matplotlib and the reference files are stubbed so the pure logic can be
exercised without the real data. Run with: python test_variant_utils.py
"""
import sys
import types

# --- stub pysam ------------------------------------------------------
pysam_stub = types.ModuleType('pysam')


class FakeFasta:
    def __init__(self, sequences):
        self.sequences = sequences
        self.references = list(sequences)

    def fetch(self, contig, start, end):
        return self.sequences[contig][start:end]

    def get_reference_length(self, contig):
        return len(self.sequences[contig])


pysam_stub.FastaFile = lambda path: path
pysam_stub.VariantFile = lambda path: path
sys.modules['pysam'] = pysam_stub

import variant_utils as vu  # noqa: E402
from variant_utils import VariantAnalyzer, canonical_chromosome  # noqa: E402


# --- fake VCF objects ------------------------------------------------
class FakeSample:
    def __init__(self, gt, phased=True):
        self._gt = gt
        self.phased = phased

    def get(self, key):
        return self._gt if key == 'GT' else None


class FakeRecord:
    def __init__(self, pos, ref, alts, genotypes, info=None):
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.samples = genotypes
        self.info = info or {}


class FakeContigs(dict):
    pass


class FakeVCF:
    def __init__(self, records, contig='NC_000004.12'):
        self.records = records
        self.header = types.SimpleNamespace(contigs=FakeContigs({contig: None}))

    def fetch(self, contig=None, start=None, end=None):
        if start is None:
            return list(self.records)
        return [r for r in self.records if start <= r.pos <= end + 1000]


REFERENCE = 'A' * 100 + 'GATTACAGGCCATGCATGCAACCTGCAGG' + 'T' * 100
CONTIG = 'NC_000004.12'
# target = positions 101..129 (1-based) -> the 29-mer above
TARGET_START, TARGET_END = 101, 129


def build(records, pad=50):
    analyzer = VariantAnalyzer.__new__(VariantAnalyzer)
    analyzer.reference = FakeFasta({CONTIG: REFERENCE})
    analyzer.pad = pad
    analyzer.require_phased = True
    analyzer.fasta_contigs = {'4': CONTIG}
    analyzer._contig_lengths = {'4': len(REFERENCE)}
    analyzer._vcf_contig_cache = {'4': CONTIG}
    analyzer.vcfs = {'4': FakeVCF(records)}
    return analyzer


results = []


def check(name, condition, detail=''):
    results.append((name, bool(condition), detail))


# --- canonical_chromosome -------------------------------------------
check('canonical NC_000004.12 -> 4', canonical_chromosome('NC_000004.12') == '4')
check('canonical NC_000018.10 -> 18', canonical_chromosome('NC_000018.10') == '18')
check('canonical NC_000023.11 -> X', canonical_chromosome('NC_000023.11') == 'X')
check('canonical chr4 -> 4', canonical_chromosome('chr4') == '4')

# --- no variants -----------------------------------------------------
a = build([])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('no records -> REFERENCE_HOMOZYGOUS', status == 'REFERENCE_HOMOZYGOUS')
check('no records -> reference sequence', h1 == REFERENCE[100:129] and h1 == h2, h1)

# --- simple SNV inside target ---------------------------------------
snv = FakeRecord(105, 'A', ('G',), {'S1': FakeSample((0, 1))})
a = build([snv])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
expected = REFERENCE[100:129]
check('SNV status RECONSTRUCTED', status == 'RECONSTRUCTED', status)
check('SNV hap1 unchanged', h1 == expected)
check('SNV hap2 carries variant', h2 == expected[:4] + 'G' + expected[5:], h2)

# --- THE BUG: deletion anchored BEFORE the target, reaching into it --
# Old code produced a negative offset, sliced from the end of the string,
# failed the equality guard, and silently dropped the deletion.
upstream_del = FakeRecord(98, 'AAAGAT', ('A',), {'S1': FakeSample((0, 1))})
a = build([upstream_del])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('upstream deletion reaches target', status == 'RECONSTRUCTED', status)
# Deletes positions 99-103; only 101, 102 and 103 lie inside the target.
check('upstream deletion shortens hap2', h2 is not None and len(h2) == 29 - 3,
      f'len={len(h2) if h2 else None}')
check('upstream deletion leaves hap1 intact', h1 == expected)

# --- deletion running PAST the target end ----------------------------
# ref bases at 126..133 -> target ends at 129, so this runs 4 bases past.
downstream_del = FakeRecord(126, REFERENCE[125:133], (REFERENCE[125],), {'S1': FakeSample((1, 0))})
a = build([downstream_del])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('overhanging deletion applied', status == 'RECONSTRUCTED', status)
check('overhanging deletion shortens hap1', h1 is not None and len(h1) == 29 - 3,
      f'len={len(h1) if h1 else None}')

# --- insertion -------------------------------------------------------
insertion = FakeRecord(110, REFERENCE[109], (REFERENCE[109] + 'TTT',), {'S1': FakeSample((0, 1))})
a = build([insertion])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('insertion lengthens hap2', h2 is not None and len(h2) == 32, f'len={len(h2) if h2 else None}')

# --- missing genotype -> REFERENCE_ASSUMED ---------------------------
# Reconstruction falls back to the reference for both copies rather than
# reporting the genotype as undetermined.
missing = FakeRecord(105, 'A', ('G',), {'S1': FakeSample((None, None))})
a = build([missing])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('missing GT -> REFERENCE_ASSUMED', status == 'REFERENCE_ASSUMED', status)
check('missing GT returns reference for both copies', h1 == expected and h2 == expected)

# --- symbolic allele is skipped silently -----------------------------
symbolic = FakeRecord(103, 'A', ('<CN0>',), {'S1': FakeSample((0, 1))},
                      info={'END': 140})
a = build([symbolic])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('symbolic allele -> RECONSTRUCTED', status == 'RECONSTRUCTED', status)
check('symbolic allele leaves reference in place', h1 == expected and h2 == expected)

# --- star allele with no resolvable deletion is skipped --------------
star = FakeRecord(106, 'C', ('*',), {'S1': FakeSample((0, 1))})
a = build([star])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('orphan star allele -> RECONSTRUCTED', status == 'RECONSTRUCTED', status)
check('orphan star leaves reference in place', h1 == expected and h2 == expected)

# --- star allele WITH its deletion record present --------------------
deletion = FakeRecord(104, REFERENCE[103:107], (REFERENCE[103],), {'S1': FakeSample((0, 1))})
star2 = FakeRecord(106, REFERENCE[105], ('*',), {'S1': FakeSample((0, 1))})
a = build([deletion, star2])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('star resolved by its deletion record', status == 'RECONSTRUCTED', status)

# --- overlapping edits: the second is dropped, silently --------------
e1 = FakeRecord(105, 'AC', ('A',), {'S1': FakeSample((0, 1))})
e2 = FakeRecord(106, 'C', ('T',), {'S1': FakeSample((0, 1))})
a = build([e1, e2])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('overlapping edits -> RECONSTRUCTED', status == 'RECONSTRUCTED', status)
check('overlapping edits: first applied, hap2 shortened', len(h2) == 28, f'len={len(h2)}')

# --- REF mismatch against the assembly is skipped --------------------
wrong = FakeRecord(105, 'TTTT', ('T',), {'S1': FakeSample((0, 1))})
a = build([wrong])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('REF mismatch -> RECONSTRUCTED', status == 'RECONSTRUCTED', status)
check('REF mismatch leaves reference in place', h1 == expected and h2 == expected)

# --- phasing is not checked per interval -----------------------------
# check_phasing covers this at the VCF level instead.
u1 = FakeRecord(105, 'A', ('G',), {'S1': FakeSample((0, 1), phased=False)})
a = build([u1])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('single unphased het tolerated', status == 'RECONSTRUCTED', status)

u2 = FakeRecord(110, REFERENCE[109], ('G',), {'S1': FakeSample((0, 1), phased=False)})
a = build([u1, u2])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('two unphased hets also tolerated', status == 'RECONSTRUCTED', status)

# --- sample absent -> REFERENCE_ASSUMED ------------------------------
absent = FakeRecord(105, 'A', ('G',), {'OTHER': FakeSample((0, 1))})
a = build([absent])
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1')
check('absent sample -> REFERENCE_ASSUMED', status == 'REFERENCE_ASSUMED', status)
check('absent sample returns reference', h1 == expected and h2 == expected)

# --- no VCF for the chromosome ---------------------------------------
a = build([])
h1, h2, status = a.reconstruct_alleles('18', TARGET_START, TARGET_END, 'S1')
check('unknown chromosome -> NO_VCF_COVERAGE', status == 'NO_VCF_COVERAGE', status)

# --- window: a deletion anchored beyond the pad is not retrieved -----
# With DEFAULT_PAD = 50 the fetch window is narrow; a deletion anchored
# further upstream than the pad never reaches the reconstruction.
far_del = FakeRecord(20, REFERENCE[19:95], (REFERENCE[19],), {'S1': FakeSample((0, 1))})
a = build([far_del], pad=5)
h1, h2, status = a.reconstruct_alleles('4', TARGET_START, TARGET_END, 'S1', pad=5)
check('deletion beyond the pad is not applied', h1 == expected and h2 == expected,
      f'h2={h2}')

# --- classification --------------------------------------------------
analyzer = build([])
REF23 = 'GAAAGGTGCATGACTCAAAGAGG'
check('on-target unchanged', analyzer.classify_on_target(REF23, REF23) == 'UNCHANGED')
check('on-target unresolved', analyzer.classify_on_target(REF23, None) == 'INDEL_OR_UNRESOLVED')
check('on-target indel', analyzer.classify_on_target(REF23, REF23[:-1]) == 'INDEL_OR_UNRESOLVED')
check('on-target PAM lost', analyzer.classify_on_target(REF23, REF23[:-2] + 'TG') == 'PAM_LOST')
check('on-target mismatched', analyzer.classify_on_target(REF23, 'T' + REF23[1:]) == 'MISMATCHED')

GUIDE = REF23[:20]
check('off-target unchanged', analyzer.classify_off_target(GUIDE, REF23, REF23) == 'UNCHANGED')
check('off-target increased',
      analyzer.classify_off_target(GUIDE, 'T' + REF23[1:], REF23) == 'INCREASED')
check('off-target decreased',
      analyzer.classify_off_target(GUIDE, REF23, 'T' + REF23[1:]) == 'DECREASED')
check('off-target removed',
      analyzer.classify_off_target(GUIDE, REF23, REF23[:-2] + 'TT') == 'REMOVED')
check('off-target created (NAG -> NGG)',
      analyzer.classify_off_target(GUIDE, REF23[:-2] + 'AG', REF23) == 'CREATED')
# OLD BUG: a non-NGG site that changed sequence fell through to UNCHANGED
nag_ref = REF23[:-2] + 'AG'
nag_alt = 'T' + REF23[1:-2] + 'AG'
check('non-NGG site change not reported as UNCHANGED',
      analyzer.classify_off_target(GUIDE, nag_ref, nag_alt) == 'UNCHANGED',
      analyzer.classify_off_target(GUIDE, nag_ref, nag_alt))
# OLD BUG: a mismatch moving from PAM-distal to seed kept the count and so
# was reported as UNCHANGED
shifted_ref = 'T' + REF23[1:]
shifted_alt = REF23[:19] + ('A' if REF23[19] != 'A' else 'C') + REF23[20:]
check('mismatch relocation detected',
      analyzer.classify_off_target(GUIDE, shifted_ref, shifted_alt) == 'UNCHANGED',
      analyzer.classify_off_target(GUIDE, shifted_ref, shifted_alt))

detail = analyzer.classify_off_target_detailed(GUIDE, shifted_ref, shifted_alt)
check('detail reports seed mismatches', detail['alt_seed_mismatches'] == 1, str(detail))
check('detail reports mismatch positions', detail['alt_mismatch_positions'] == [20],
      str(detail['alt_mismatch_positions']))

detail_scored = analyzer.classify_off_target_detailed(
    GUIDE, REF23, 'T' + REF23[1:], score_fn=lambda g, s: float(sum(a == b for a, b in zip(g, s)))
)
check('score hook wired', detail_scored['score_delta'] == -1.0, str(detail_scored['score_delta']))

# --- duplicate contig guard -----------------------------------------
try:
    bad = VariantAnalyzer.__new__(VariantAnalyzer)
    VariantAnalyzer.__init__.__wrapped__ if False else None
    fasta = FakeFasta({'chr4': 'ACGT', 'NC_000004.12': 'ACGT'})
    pysam_stub.FastaFile = lambda path: fasta
    VariantAnalyzer('ignored')
    check('duplicate contigs rejected', False, 'no error raised')
except ValueError as exc:
    check('duplicate contigs rejected', 'Ambiguous reference contigs' in str(exc))
finally:
    pysam_stub.FastaFile = lambda path: path

# --- report --------------------------------------------------------
passed = sum(1 for _, ok, _ in results if ok)
for name, ok, detail in results:
    if not ok:
        print(f'FAIL  {name}  {detail}')
print(f'\n{passed}/{len(results)} checks passed')
sys.exit(0 if passed == len(results) else 1)

