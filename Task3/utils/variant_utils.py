import re
from pathlib import Path
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

# =====================================================================
# Sequence and coordinate helpers
# =====================================================================

def reverse_complement(sequence):
    return sequence.upper().translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def canonical_chromosome(value):
    """Normalise a contig name to a bare chromosome label ('4', '18', 'X')."""
    text = str(value).strip()
    match = re.fullmatch(r'chr([0-9]{1,2}|X|Y|M|MT)', text, flags=re.IGNORECASE)
    if match:
        token = match.group(1).upper()
        return 'M' if token in ('M', 'MT') else token
    match = re.fullmatch(r'NC_0000(\d{2})\.\d+', text)
    if match:
        number = int(match.group(1))
        if number == 23:
            return 'X'
        if number == 24:
            return 'Y'
        return str(number)
    return text


def orient_sequence(genomic_sequence, strand):
    if genomic_sequence is None:
        return None
    return genomic_sequence if strand == '+' else reverse_complement(genomic_sequence)


def record_end(record):
    return int(record.info['END']) if 'END' in record.info else record.pos + len(record.ref) - 1


# How far either side of the target interval to fetch reference sequence and
# variant records.
DEFAULT_PAD = 50


class VariantAnalyzer:
    """Reconstructs personalised haplotype sequences from a reference assembly
    and a phased VCF.

    Where a genotype is missing or a variant cannot be applied, the reference
    sequence is used for that position. Use `check_phasing` to confirm the
    input is phased before relying on haplotype-level results.
    """

    def __init__(self, fasta_path, pad=DEFAULT_PAD, require_phased=True):
        self.reference = pysam.FastaFile(str(fasta_path))
        self.vcfs = {}
        self.pad = int(pad)
        self.require_phased = bool(require_phased)

        self.fasta_contigs = {}
        collisions = {}
        for contig in self.reference.references:
            key = canonical_chromosome(contig)
            if key in self.fasta_contigs:
                collisions.setdefault(key, [self.fasta_contigs[key]]).append(contig)
            else:
                self.fasta_contigs[key] = contig
        if collisions:
            raise ValueError(
                'Ambiguous reference contigs — several FASTA sequences map to the '
                'same chromosome label, so the wrong one could be used: '
                + '; '.join(f'{k}: {v}' for k, v in collisions.items())
            )

        self._contig_lengths = {
            canonical_chromosome(c): self.reference.get_reference_length(c)
            for c in self.reference.references
            if canonical_chromosome(c) in self.fasta_contigs
        }
        self._vcf_contig_cache = {}

    # -- setup ---------------------------------------------------------

    def add_vcf(self, chrom, vcf_path):
        key = canonical_chromosome(chrom)
        vcf = pysam.VariantFile(str(vcf_path))
        self.vcfs[key] = vcf
        self._vcf_contig_cache.pop(key, None)
        return vcf

    def vcf_contig(self, chrom):
        """The contig name as spelled inside the VCF for a canonical label."""
        chrom = canonical_chromosome(chrom)
        if chrom in self._vcf_contig_cache:
            return self._vcf_contig_cache[chrom]
        vcf = self.vcfs[chrom]
        resolved = None
        for contig in vcf.header.contigs:
            if canonical_chromosome(contig) == chrom:
                resolved = contig
                break
        self._vcf_contig_cache[chrom] = resolved
        return resolved

    def check_phasing(self, chrom, max_records=5000):
        """Report what fraction of heterozygous genotypes are phased.

        Haplotype-level results are only meaningful on phased data, so call
        this once per VCF and record the answer in the write-up.
        """
        chrom = canonical_chromosome(chrom)
        vcf = self.vcfs[chrom]
        het = phased = 0
        for i, record in enumerate(vcf.fetch()):
            if i >= max_records:
                break
            for sample in record.samples.values():
                gt = sample.get('GT')
                if gt is None or len(gt) < 2 or any(a is None for a in gt[:2]):
                    continue
                if gt[0] != gt[1]:
                    het += 1
                    if sample.phased:
                        phased += 1
        return {
            'chrom': chrom,
            'records_examined': min(i + 1, max_records),
            'heterozygous_genotypes': het,
            'phased_heterozygous': phased,
            'phased_fraction': (phased / het) if het else float('nan'),
        }

    # -- record retrieval ----------------------------------------------

    def _get_overlapping_records(self, vcf, chrom, start, end, pad=None):
        """Every record whose REF span intersects the 1-based interval
        [start, end]. Fetches `pad` bases either side so that deletions
        anchored just outside the interval but reaching into it are not
        missed; a deletion anchored further out than `pad` is not retrieved."""
        pad = self.pad if pad is None else pad
        fetch_start = max(0, start - 1 - pad)
        fetch_end = end + pad
        return [
            r for r in vcf.fetch(chrom, fetch_start, fetch_end)
            if r.pos <= end and record_end(r) >= start
        ]

    # -- haplotype reconstruction --------------------------------------

    def reconstruct_alleles(self, chrom, start, end, sample_id, pad=None):
        """Return (allele_1, allele_2, status) for the 1-based interval
        [start, end].

        `status` is one of NO_VCF_COVERAGE, REFERENCE_HOMOZYGOUS,
        REFERENCE_ASSUMED or RECONSTRUCTED. A missing or unusable genotype
        falls back to the reference sequence for both copies and is reported
        as REFERENCE_ASSUMED; every other situation that cannot be applied is
        skipped silently, leaving reference sequence at that position.

        Variants are applied on a per-reference-base grid spanning a padded
        window, so deletions that start before the interval or run past its
        end shorten the returned sequence correctly instead of being dropped.
        """
        chrom = canonical_chromosome(chrom)
        if chrom not in self.vcfs:
            return None, None, 'NO_VCF_COVERAGE'

        pad = self.pad if pad is None else pad
        fasta_chrom = self.fasta_contigs.get(chrom, chrom)
        contig_length = self._contig_lengths.get(chrom)

        window_start = max(1, start - pad)
        window_end = end + pad
        if contig_length:
            window_end = min(window_end, contig_length)

        padded = self.reference.fetch(fasta_chrom, window_start - 1, window_end).upper()
        left = start - window_start
        right = end - window_start + 1
        target_reference = padded[left:right]

        records = self._get_overlapping_records(
            self.vcfs[chrom], self.vcf_contig(chrom), start, end, pad
        )
        if not records:
            return target_reference, target_reference, 'REFERENCE_HOMOZYGOUS'

        # One cell per reference base. A substitution replaces a cell, a
        # deletion empties the cells it removes, an insertion puts a longer
        # string in the anchor cell.
        cells = [list(padded), list(padded)]
        touched = [[False] * len(padded) for _ in range(2)]

        for record in records:
            try:
                sample = record.samples[sample_id]
            except (KeyError, IndexError):
                return target_reference, target_reference, 'REFERENCE_ASSUMED'

            genotype = sample.get('GT')
            if genotype is None or len(genotype) < 2 or any(a is None for a in genotype[:2]):
                return target_reference, target_reference, 'REFERENCE_ASSUMED'

            for haplotype, allele_index in enumerate(genotype[:2]):
                if allele_index == 0:
                    continue
                if record.alts is None or allele_index > len(record.alts):
                    continue

                alt = str(record.alts[allele_index - 1]).upper()
                reference_allele = str(record.ref).upper()
                offset = record.pos - window_start

                # Symbolic alternates, spanning-deletion '*' alleles, edits
                # falling outside the window, reference mismatches and edits
                # overlapping one already applied are all skipped, leaving
                # reference sequence in place.
                if alt.startswith('<') or alt == '*':
                    continue
                if offset < 0 or offset + len(reference_allele) > len(padded):
                    continue
                if padded[offset:offset + len(reference_allele)] != reference_allele:
                    continue
                if any(touched[haplotype][offset:offset + len(reference_allele)]):
                    continue

                cells[haplotype][offset] = alt
                for k in range(1, len(reference_allele)):
                    cells[haplotype][offset + k] = ''
                for k in range(len(reference_allele)):
                    touched[haplotype][offset + k] = True

        alleles = [''.join(cells[haplotype][left:right]) for haplotype in range(2)]
        return alleles[0], alleles[1], 'RECONSTRUCTED'

    # -- classification ------------------------------------------------

    @staticmethod
    def _pam_class(sequence):
        if sequence is None or len(sequence) < 3:
            return 'NONE'
        if sequence[-2:] == 'GG':
            return 'NGG'
        if sequence[-2:] == 'AG':
            return 'NAG'
        return 'NONE'

    def classify_on_target(self, ref_23mer, alt_23mer):
        """INDEL_OR_UNRESOLVED | PAM_LOST | MISMATCHED | UNCHANGED."""
        if not alt_23mer or len(alt_23mer) != len(ref_23mer):
            return 'INDEL_OR_UNRESOLVED'
        if self._pam_class(alt_23mer) != 'NGG':
            return 'PAM_LOST'
        protospacer_length = len(ref_23mer) - 3
        mismatches = sum(
            a != b for a, b in zip(ref_23mer[:protospacer_length], alt_23mer[:protospacer_length])
        )
        return 'UNCHANGED' if mismatches == 0 else 'MISMATCHED'

    def classify_off_target_detailed(self, guide_20mer, ref_site, alt_site, score_fn=None):
        """Full comparison of a candidate off-target site between the
        reference and one personalised haplotype.

        Returns a dict with the effect label plus the mismatch counts and
        positions behind it, so downstream analysis can weight a seed-region
        change differently from a PAM-distal one.

        `score_fn(guide_20mer, site_23mer) -> float` is an optional hook for a
        position-weighted model such as CFD; when supplied, the scores and
        their delta are included.
        """
        result = {
            'effect': None,
            'ref_pam': self._pam_class(ref_site),
            'alt_pam': self._pam_class(alt_site),
            'ref_mismatches': None,
            'alt_mismatches': None,
            'alt_mismatch_positions': None,
            'alt_seed_mismatches': None,
            'ref_score': None,
            'alt_score': None,
            'score_delta': None,
        }

        if alt_site is None or len(alt_site) != len(ref_site):
            result['effect'] = 'UNRESOLVED'
            return result

        protospacer_length = len(guide_20mer)
        ref_mismatch_positions = [
            i + 1 for i, (a, b) in enumerate(zip(guide_20mer, ref_site[:protospacer_length])) if a != b
        ]
        alt_mismatch_positions = [
            i + 1 for i, (a, b) in enumerate(zip(guide_20mer, alt_site[:protospacer_length])) if a != b
        ]
        result['ref_mismatches'] = len(ref_mismatch_positions)
        result['alt_mismatches'] = len(alt_mismatch_positions)
        result['alt_mismatch_positions'] = alt_mismatch_positions
        # Seed = the 10 nucleotides proximal to the PAM, where mismatches are
        # least tolerated.
        result['alt_seed_mismatches'] = sum(1 for p in alt_mismatch_positions if p > protospacer_length - 10)

        if score_fn is not None:
            result['ref_score'] = score_fn(guide_20mer, ref_site)
            result['alt_score'] = score_fn(guide_20mer, alt_site)
            result['score_delta'] = result['alt_score'] - result['ref_score']

        ref_ngg = result['ref_pam'] == 'NGG'
        alt_ngg = result['alt_pam'] == 'NGG'

        if not ref_ngg and alt_ngg:
            result['effect'] = 'CREATED'
        elif ref_ngg and not alt_ngg:
            result['effect'] = 'REMOVED'
        elif alt_ngg and result['alt_mismatches'] < result['ref_mismatches']:
            result['effect'] = 'INCREASED'
        elif alt_ngg and result['alt_mismatches'] > result['ref_mismatches']:
            result['effect'] = 'DECREASED'
        else:
            # Covers an unchanged site, a site with no canonical PAM either
            # way, and one whose mismatch count is unchanged.
            result['effect'] = 'UNCHANGED'
        return result

    def classify_off_target(self, guide_20mer, ref_site, alt_site, score_fn=None):
        return self.classify_off_target_detailed(guide_20mer, ref_site, alt_site, score_fn)['effect']


# =====================================================================
# Reporting helpers
# =====================================================================

def targetability_counts(frame, status_column='target_status'):
    """Per-allele counts split into targetable / disrupted / unknown.

    Keeping 'unknown' separate from 'disrupted' stops missing data from being
    reported as a biological finding.
    """
    categories = {
        'UNCHANGED': 'targetable',
        'MISMATCHED': 'disrupted',
        'PAM_LOST': 'disrupted',
        'INDEL_OR_UNRESOLVED': 'disrupted',
    }
    bucket = frame[status_column].map(categories).fillna('unknown')
    return (
        frame.assign(bucket=bucket)
        .groupby(['gene', 'guide_id', 'bucket'])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )


# =====================================================================
# Genome-wide off-target visualisation
# =====================================================================

# Colours are taken from the palettes already used elsewhere in these
# notebooks, so the genome-wide figures sit alongside the on-target heatmaps
# and the off-target risk breakdown without a style change.
#
# Red/green was the obvious pairing from the on-target heatmap, but it fails
# colour-vision separation outright (OKLab dE 1.2 under deuteranopia - the two
# bars are effectively one colour). The red is kept for risk increase and
# paired with the blue already used for REMOVED in the off-target breakdown:
# that pair passes every check (CVD dE 14.9 protan, normal-vision 29.7,
# both above 3:1 on white).
RISK_UP_COLOUR = '#E45756'      # red, as in the on-target availability heatmap
RISK_DOWN_COLOUR = '#1f77b4'    # blue, as in the off-target risk breakdown
ANCESTRY_CMAP = 'YlOrRd'        # as in the ancestry disruption heatmap

RISK_UP_EFFECTS = ('CREATED', 'INCREASED')
RISK_DOWN_EFFECTS = ('REMOVED', 'DECREASED')
NEUTRAL_EFFECTS = ('UNRESOLVED',)


def chromosome_sort_key(value):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def _render_events_by_chromosome(frame, chroms, output_path, heading, show):
    """Overview-and-detail figure of off-target events per chromosome.

    The upper panel is the honest overview: both directions on one shared,
    zero-based linear scale, so the relative magnitudes are exactly as the data
    has them. Risk-increasing counts are one to two orders of magnitude smaller
    and are hairlines at that scale, so the lower panel magnifies them.

    This works where two side-by-side panels on independent scales did not: the
    lower panel is an explicit zoom on a sub-range of the same axis, labelled as
    such, rather than a second category whose bar lengths invite comparison
    with the first. A log axis was rejected - bar length encodes magnitude
    proportionally, so a log scale misstates ratios, and several chromosomes
    have a count of zero.
    """
    positions = np.arange(len(chroms))
    up = (frame[frame['offtarget_effect'].isin(RISK_UP_EFFECTS)]
          .groupby('chrom').size().reindex(chroms, fill_value=0))
    down = (frame[frame['offtarget_effect'].isin(RISK_DOWN_EFFECTS)]
            .groupby('chrom').size().reindex(chroms, fill_value=0))

    fig, (ax_all, ax_zoom) = plt.subplots(
        2, 1, figsize=(12, 9), dpi=300, sharex=True,
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.28})
    chrom_labels = [f'chr{c}' for c in chroms]

    # -- overview: both directions, one shared scale --------------------
    ax_all.bar(positions, up.to_numpy(), width=0.62, color=RISK_UP_COLOUR,
               label='Risk increased (created / increased)')
    ax_all.bar(positions, -down.to_numpy(), width=0.62, color=RISK_DOWN_COLOUR,
               label='Risk decreased (removed / decreased)')
    ax_all.axhline(0, color='black', linewidth=0.8)
    span = max(up.max(), down.max()) or 1
    ax_all.set_ylim(-span * 1.15, span * 1.15)
    # Only the risk-decreasing bars are annotated here. The risk-increasing
    # counts are hairlines at this scale, so their labels would sit on the zero
    # line and collide with these; the detail panel below already carries them.
    for x, v in zip(positions, down.to_numpy()):
        if v:
            ax_all.annotate(f'{v:,}', (x, -v), xytext=(0, -4), textcoords='offset points',
                            ha='center', va='top', fontsize=7.5, weight='bold')
    ax_all.set_ylabel('Observed Events\n(Haplotype Count)', fontsize=10, weight='bold')
    ax_all.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{abs(int(v)):,}'))
    ax_all.yaxis.grid(True, linestyle=':', linewidth=0.6, alpha=0.6)
    ax_all.set_axisbelow(True)
    # 'best' rather than a fixed corner: which corner is free depends on the
    # gene. A fixed 'lower right' buried SERPINB2's chr16 and chr17 bars and
    # their labels under the legend box.
    ax_all.legend(loc='best', fontsize=9, frameon=True, framealpha=0.95)
    ax_all.set_title('Panel A: Both directions, shared scale',
                     fontsize=10, weight='bold', pad=6)
    # The risk-increasing bars are unreadable at this scale, so the panel says
    # where to read them rather than leaving the reader to infer it. Grey and
    # unemphasised: it is a pointer, not a finding.
    ax_all.text(0.015, 0.97, 'See Panel B for risk-increasing event counts',
                transform=ax_all.transAxes, ha='left', va='top',
                fontsize=8.5, style='italic', color='#666666')
    # The axes are shared, so matplotlib would label only the lower one. The
    # bars here hang from a zero line in the middle of the panel, well away
    # from the axis, so without its own labels this panel has to be read by
    # counting bars across from Panel B.
    ax_all.set_xticks(positions)
    ax_all.set_xticklabels(chrom_labels, rotation=45, ha='right',
                           fontsize=8, weight='bold')
    ax_all.tick_params(axis='x', labelbottom=True)

    # -- detail: the risk-increasing bars, magnified ---------------------
    ax_zoom.bar(positions, up.to_numpy(), width=0.62, color=RISK_UP_COLOUR)
    zoom_span = up.max() or 1
    ax_zoom.set_ylim(0, zoom_span * 1.22)
    for x, v in zip(positions, up.to_numpy()):
        if v:
            ax_zoom.annotate(f'{v:,}', (x, v), xytext=(0, 3), textcoords='offset points',
                             ha='center', fontsize=8, weight='bold')
    ax_zoom.set_ylabel('Risk-Increasing\nEvents', fontsize=10, weight='bold')
    ax_zoom.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax_zoom.yaxis.grid(True, linestyle=':', linewidth=0.6, alpha=0.6)
    ax_zoom.set_axisbelow(True)
    # Where the two directions are already comparable in size the lower panel
    # is not a magnification of anything - SERPINB2 comes out at 1x - so the
    # claim is only made when there is a real difference in scale.
    magnification = span / zoom_span
    zoom_title = 'Panel B: Risk-increasing events only'
    if magnification >= 1.5:
        zoom_title = f'{zoom_title}, magnified ({magnification:.0f}x the scale in Panel A)'
    ax_zoom.set_title(zoom_title, fontsize=10, weight='bold', pad=6)
    ax_zoom.set_xticks(positions)
    ax_zoom.set_xticklabels(chrom_labels, rotation=45, ha='right',
                            fontsize=9, weight='bold')
    ax_zoom.set_xlabel('Chromosome', fontsize=11, weight='bold')

    fig.suptitle(heading, fontsize=13, weight='bold')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if output_path:
        fig.savefig(output_path, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)
    return pd.DataFrame({'risk_increased': up, 'risk_decreased': down})


def plot_offtarget_events_by_chromosome(events, output_path=None, title=None,
                                        show=True, split_by_gene=True):
    """Diverging horizontal bar chart of off-target events per chromosome.

    Risk-increasing events (CREATED, INCREASED) extend right of the zero line
    and risk-decreasing events (REMOVED, DECREASED) extend left, so direction
    is carried by position as well as by colour. Counts are haplotype
    observations.

    With `split_by_gene` and a 'gene' column covering more than one gene, each
    gene is written as its own figure, with the gene appended to the file name
    and its own count scale, so the two can be reported separately.

    `events` needs 'chrom' and 'offtarget_effect' columns and may contain any
    effect label; UNCHANGED rows are ignored. Non-directional effects
    (UNRESOLVED) have no direction to plot, so they are counted in the
    returned frame's source data but carry no bar of their own.
    """
    frame = events.copy()
    frame['chrom'] = frame['chrom'].astype(str)
    cohort_size = frame['sample_id'].nunique() if 'sample_id' in frame.columns else None
    chroms = sorted(frame['chrom'].unique(), key=chromosome_sort_key)

    genes = (sorted(frame['gene'].unique())
             if split_by_gene and 'gene' in frame.columns and frame['gene'].nunique() > 1
             else [None])

    results = {}
    for gene in genes:
        subset = frame if gene is None else frame[frame['gene'] == gene]
        heading = title or 'Count of Predicted Off-Target Risk Events by Chromosome'
        if gene:
            heading = f'{gene}: {heading}'
        if cohort_size and 'N =' not in heading:
            heading = f'{heading}\n(N = {cohort_size:,} Individuals)'

        path = output_path
        if path is not None and gene:
            path = Path(path)
            path = path.with_name(f'{path.stem}_{gene.lower()}{path.suffix}')

        results[gene or 'ALL'] = _render_events_by_chromosome(
            subset, chroms, path, heading, show)

    return pd.concat(results, names=['gene', 'chrom'])


def risk_increasing_site_report(events, validated_offtargets=None, cohort_size=None,
                                metadata=None, population_order=('EUR', 'AFR', 'EAS', 'SAS', 'AMR')):
    """One row per candidate site where some individual's personal sequence
    raises predicted off-target risk.

    This is the table behind the by-guide and by-chromosome figures: it names
    the guide, the chromosome and position, the locus, CRISPOR's reference
    scores, and how many individuals carry the risk-raising allele. Pass
    `metadata` (columns 'sample_id' and 'super_population') to add a carrier
    percentage per super-population.

    Returns the per-site table sorted by carrier count, descending.
    """
    risk = events[events['offtarget_effect'].isin(RISK_UP_EFFECTS)].copy()
    if risk.empty:
        return pd.DataFrame()
    risk['chrom'] = risk['chrom'].astype(str)

    keys = ['gene', 'guide_id', 'offtarget_id', 'chrom', 'start_1based']
    report = (risk.groupby(keys)
              .agg(observations=('offtarget_effect', 'size'),
                   carriers=('sample_id', 'nunique'))
              .reset_index())

    if cohort_size:
        report['carrier_pct'] = (100 * report['carriers'] / cohort_size).round(3)

    if validated_offtargets is not None:
        extra = [c for c in ['strand', 'locus_description', 'reference_mismatch_count_crispor',
                             'reference_mit_score', 'reference_cfd_score']
                 if c in validated_offtargets.columns]
        report = report.merge(validated_offtargets[['offtarget_id'] + extra],
                              on='offtarget_id', how='left')

    if metadata is not None and 'super_population' in metadata.columns:
        sizes = metadata['super_population'].value_counts().to_dict()
        pops = [p for p in population_order if sizes.get(p)]
        labelled = risk.merge(metadata[['sample_id', 'super_population']].drop_duplicates(),
                              on='sample_id', how='left', suffixes=('', '_meta'))
        column = ('super_population_meta' if 'super_population_meta' in labelled.columns
                  else 'super_population')
        per_pop = (labelled.drop_duplicates(['offtarget_id', 'sample_id'])
                   .groupby(['offtarget_id', column]).size().unstack(fill_value=0)
                   .reindex(columns=pops, fill_value=0))
        per_pop = (per_pop.divide([sizes[p] for p in pops], axis=1) * 100).round(3)
        per_pop.columns = [f'{p}_carrier_pct' for p in pops]
        report = report.merge(per_pop.reset_index(), on='offtarget_id', how='left')

    return report.sort_values('carriers', ascending=False).reset_index(drop=True)


def risk_increasing_guide_by_chromosome(events):
    """Matrix of risk-increasing observations: one row per gene/guide, one
    column per chromosome. Shows at a glance which chromosomes a guide's
    off-target burden comes from."""
    risk = events[events['offtarget_effect'].isin(RISK_UP_EFFECTS)].copy()
    if risk.empty:
        return pd.DataFrame()
    risk['chrom'] = risk['chrom'].astype(str)
    matrix = (risk.groupby(['gene', 'guide_id', 'chrom']).size()
              .unstack(fill_value=0))
    matrix = matrix[sorted(matrix.columns, key=chromosome_sort_key)]
    matrix.columns = [f'chr{c}' for c in matrix.columns]
    return matrix.assign(total=matrix.sum(axis=1)).sort_values('total', ascending=False)


def plot_offtarget_events_by_guide(events, output_path=None, title=None,
                                   show=True, split_by_gene=True,
                                   all_guides=None):
    """Diverging horizontal bar chart of off-target events per guide.

    Same encoding as the by-chromosome figure - risk-increasing right of zero,
    risk-decreasing left - but grouped by guide, so the guides responsible for
    a gene's off-target burden can be read off directly. Guides are ordered by
    risk-increasing count, and a guide with none is shown at zero rather than
    dropped, since "this guide is clean" is itself a result.

    The number of distinct candidate sites behind each risk-increasing bar is
    annotated in brackets: many observations at one site is a common variant,
    whereas the same count spread over many sites is a different problem.

    `all_guides` is the guide manifest (columns 'gene' and 'guide_id'). Pass it
    whenever `events` holds only altered observations, as the N=2,504 notebook's
    table does - otherwise a guide with no events is missing from the figure
    entirely, and a clean guide would be indistinguishable from an absent one.
    """
    frame = events.copy()
    cohort_size = frame['sample_id'].nunique() if 'sample_id' in frame.columns else None
    genes = (sorted(frame['gene'].unique())
             if split_by_gene and 'gene' in frame.columns and frame['gene'].nunique() > 1
             else [None])

    results = {}
    for gene in genes:
        subset = frame if gene is None else frame[frame['gene'] == gene]
        if all_guides is not None:
            catalogue = all_guides if gene is None else all_guides[all_guides['gene'] == gene]
            guides = sorted(catalogue['guide_id'].unique())
        else:
            guides = sorted(subset['guide_id'].unique())
        up_rows = subset[subset['offtarget_effect'].isin(RISK_UP_EFFECTS)]
        down_rows = subset[subset['offtarget_effect'].isin(RISK_DOWN_EFFECTS)]
        up = up_rows.groupby('guide_id').size().reindex(guides, fill_value=0)
        down = down_rows.groupby('guide_id').size().reindex(guides, fill_value=0)
        sites = up_rows.groupby('guide_id')['offtarget_id'].nunique().reindex(guides, fill_value=0)
        order = up.sort_values(ascending=True).index
        up, down, sites = up[order], down[order], sites[order]

        positions = np.arange(len(order))
        plt.figure(figsize=(10, max(3.5, 0.55 * len(order) + 2.5)), dpi=300)
        ax = plt.gca()
        ax.barh(positions, up.to_numpy(), height=0.62, color=RISK_UP_COLOUR,
                label='Risk increased (created / increased)')
        ax.barh(positions, -down.to_numpy(), height=0.62, color=RISK_DOWN_COLOUR,
                label='Risk decreased (removed / decreased)')

        span = max(up.max(), down.max()) or 1
        for y, rise, fall, n_sites in zip(positions, up.to_numpy(), down.to_numpy(), sites.to_numpy()):
            label = f'{rise:,} ({n_sites} site{"s" if n_sites != 1 else ""})' if rise else '0'
            ax.annotate(label, (rise, y), xytext=(4, 0), textcoords='offset points',
                        ha='left', va='center', fontsize=9, weight='bold')
            if fall:
                ax.annotate(f'{fall:,}', (-fall, y), xytext=(-4, 0), textcoords='offset points',
                            ha='right', va='center', fontsize=9, weight='bold')

        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_yticks(positions)
        ax.set_yticklabels(order, fontsize=10, weight='bold')
        ax.set_xlim(-span * 1.30, span * 1.30)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{abs(int(v)):,}'))
        ax.xaxis.grid(True, linestyle=':', linewidth=0.6, alpha=0.6)
        ax.set_axisbelow(True)

        heading = title or 'Off-target Risk Count by Guide'
        if gene:
            heading = f'{gene}: {heading}'
        if cohort_size and 'N =' not in heading:
            heading = f'{heading}\n(N = {cohort_size:,} Individuals)'
        plt.title(heading, fontsize=13, weight='bold', pad=15)
        plt.xlabel('Observed Events (Haplotype Count)', fontsize=11, weight='bold')
        plt.ylabel('gRNA Identifier', fontsize=11, weight='bold')
        ax.legend(loc='lower right', fontsize=9, frameon=True)
        plt.tight_layout()

        path = output_path
        if path is not None and gene:
            path = Path(path)
            path = path.with_name(f'{path.stem}_{gene.lower()}{path.suffix}')
        if path:
            plt.savefig(path, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()

        results[gene or 'ALL'] = pd.DataFrame(
            {'risk_increased': up, 'risk_decreased': down, 'risk_increasing_sites': sites})

    return pd.concat(results, names=['gene', 'guide_id'])


def plot_risk_site_ancestry_heatmap(events, population_sizes, output_path=None,
                                    title=None, show=True, top_n=15,
                                    population_order=('EUR', 'AFR', 'EAS', 'SAS', 'AMR')):
    """Heatmap of carrier percentage per super-population at each
    risk-increasing off-target site.

    A carrier is an individual with at least one haplotype whose personal
    sequence raises predicted risk at that site. The denominator is the whole
    super-population, so a zero cell means nobody in that group carries it.
    Every cell is annotated, so magnitude is never carried by colour alone.

    `events` needs 'super_population', 'sample_id', 'gene', 'guide_id',
    'chrom', 'start_1based' and 'offtarget_effect'; `population_sizes` maps
    super-population code to its cohort size.
    """
    risk = events[events['offtarget_effect'].isin(RISK_UP_EFFECTS)].copy()
    if risk.empty:
        return None, pd.DataFrame()

    risk['site'] = (risk['gene'] + ' ' + risk['guide_id']
                    + ' | chr' + risk['chrom'].astype(str) + ':' + risk['start_1based'].astype(str))
    # A population with no individuals has no denominator, so it is left out
    # rather than producing a NaN cell.
    pops = [p for p in population_order if population_sizes.get(p)]
    # Two lines per label: with five super-populations the single-line form
    # overruns its column and the names collide.
    labels = {'EUR': 'European\n(EUR)', 'AFR': 'African\n(AFR)', 'EAS': 'East Asian\n(EAS)',
              'SAS': 'South Asian\n(SAS)', 'AMR': 'Admixed Amer.\n(AMR)'}

    carriers = (risk.drop_duplicates(['site', 'sample_id', 'super_population'])
                .groupby(['site', 'super_population']).size().unstack(fill_value=0)
                .reindex(columns=pops, fill_value=0))
    percent = carriers.divide([population_sizes[p] for p in pops], axis=1) * 100
    percent = percent.loc[percent.max(axis=1).sort_values(ascending=False).index]
    # Both frames are returned together, so they are kept in the same row
    # order: the caller writes them as a matched pair of tables.
    carriers = carriers.reindex(index=percent.index)

    # At population scale most sites are carried by a handful of people and
    # their rows are a wall of zeros, so only the most frequent are plotted.
    # The truncation applies to the figure alone; the returned tables keep
    # every site, so the full table is written alongside the figure.
    total_sites = len(percent)
    truncated = top_n is not None and total_sites > top_n
    plotted = percent.head(top_n) if truncated else percent

    plt.figure(figsize=(10.5, max(3.0, 0.75 * len(plotted) + 1.8)), dpi=300)
    ax = sns.heatmap(
        plotted.rename(columns=labels), annot=True, fmt='.2f', cmap=ANCESTRY_CMAP, vmin=0,
        cbar_kws={'label': 'Carriers of \u2265 1 Risk-Increasing Allele (%)'},
        linewidths=0.5, linecolor='lightgray', annot_kws={'weight': 'bold'},
    )
    # Same convention as the by-chromosome figure and the other report
    # figures: the cohort size is stated in the title and derived from the
    # data, so it cannot drift from the cohort actually plotted.
    cohort_size = sum(population_sizes[p] for p in pops)
    heading = title or 'Risk-Increasing Off-Target Sites by Genetic Ancestry'
    if truncated:
        heading = f'{heading}\nTop {len(plotted)} of {total_sites} Sites by Carrier Frequency'
    if cohort_size and 'N =' not in heading:
        heading = f'{heading}\n(N = {cohort_size:,} Individuals)'
    plt.title(heading, fontsize=13, weight='bold', pad=15)
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(fontsize=10, weight='bold')
    plt.yticks(rotation=0, fontsize=9)
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()
    return percent, carriers


# =====================================================================
# Genomic Architecture and gRNA Target Visualisation
# =====================================================================

GENE_MODELS = {
    'CXCL11': {
        'gene_name': 'CXCL11',
        'chrom': '4',
        'chrom_accession': 'NC_000004.12',
        'transcript': 'NM_005409.5',
        'strand': '-',
        'gene_start': 76033682,
        'gene_end': 76041415,
        'mode': 'CRISPRi',
        'exons': [
            {'exon_num': 1, 'start': 76035927, 'end': 76036197, 'is_first': True},
            {'exon_num': 2, 'start': 76035216, 'end': 76035342, 'is_first': False},
            {'exon_num': 3, 'start': 76035047, 'end': 76035119, 'is_first': False},
            {'exon_num': 4, 'start': 76033682, 'end': 76034816, 'is_first': False},
        ]
    },
    'SERPINB2': {
        'gene_name': 'SERPINB2',
        'chrom': '18',
        'chrom_accession': 'NC_000018.10',
        'transcript': 'NM_001143818.1',
        'strand': '+',
        'gene_start': 63871692,
        'gene_end': 63903888,
        'mode': 'CRISPRa',
        'exons': [
            {'exon_num': 1, 'start': 63887705, 'end': 63887770, 'is_first': True},
            {'exon_num': 2, 'start': 63889838, 'end': 63890095, 'is_first': False},
            {'exon_num': 3, 'start': 63891436, 'end': 63891612, 'is_first': False},
            {'exon_num': 4, 'start': 63895264, 'end': 63895383, 'is_first': False},
            {'exon_num': 5, 'start': 63897091, 'end': 63897219, 'is_first': False},
            {'exon_num': 6, 'start': 63897727, 'end': 63897844, 'is_first': False},
            {'exon_num': 7, 'start': 63901740, 'end': 63901882, 'is_first': False},
            {'exon_num': 8, 'start': 63902404, 'end': 63902568, 'is_first': False},
            {'exon_num': 9, 'start': 63902901, 'end': 63903890, 'is_first': False},
        ]
    }
}


def plot_gene_grna_architecture(
    gene_or_chrom,
    guides_df=None,
    figsize=(15, 10.5),
    output_path=None,
    show=True
):
    """
    Visualises full gene sequence architecture, highlights the first exon, and plots
    selected gRNA protospacers with their respective PAM regions with zero overlap.

    Parameters
    ----------
    gene_or_chrom : str or int
        Gene symbol ('CXCL11', 'SERPINB2') or chromosome ('4', '18').
    guides_df : pd.DataFrame, optional
        DataFrame containing guide metadata. If None, dynamically loaded from manifest.
    figsize : tuple, default (15, 10.5)
        Figure dimensions (width, height).
    output_path : str or Path, optional
        Filepath where the rendered figure will be saved.
    show : bool, default True
        Whether to display the figure with plt.show().

    Returns
    -------
    matplotlib.figure.Figure
        The rendered matplotlib Figure.
    """
    key = str(gene_or_chrom).strip().upper()
    if key in ('4', 'CHR4', 'CXCL11'):
        gene = 'CXCL11'
    elif key in ('18', 'CHR18', 'SERPINB2'):
        gene = 'SERPINB2'
    else:
        raise ValueError(f"Unsupported gene/chrom '{gene_or_chrom}'. Expected 'CXCL11'/'4' or 'SERPINB2'/'18'.")

    model = GENE_MODELS[gene]

    # Attempt to load guides if not supplied
    if guides_df is None:
        for candidate_path in [
            Path('results/tables/task3_guide_manifest.tsv'),
            Path('../results/tables/task3_guide_manifest.tsv'),
            Path(f'results/tables/task2_{gene.lower()}_{model["mode"].lower()}_final_pool.csv'),
            Path(f'../results/tables/task2_{gene.lower()}_{model["mode"].lower()}_final_pool.csv'),
        ]:
            if candidate_path.exists():
                sep = '\t' if candidate_path.suffix in ('.tsv', '.txt') else ','
                guides_df = pd.read_csv(candidate_path, sep=sep)
                break

    if guides_df is None:
        raise FileNotFoundError(f"No guide data supplied and manifest not found for {gene}.")

    # Filter guides for current target gene/chromosome
    if 'gene' in guides_df.columns:
        gene_guides = guides_df[guides_df['gene'] == gene].copy()
    elif 'fasta_chrom' in guides_df.columns:
        gene_guides = guides_df[guides_df['fasta_chrom'].astype(str) == model['chrom']].copy()
    else:
        gene_guides = guides_df.copy()

    # Identify relevant columns
    id_col = 'guide_id' if 'guide_id' in gene_guides.columns else '#guideId'
    strand_col = 'genomic_strand' if 'genomic_strand' in gene_guides.columns else 'mapped_orientation'

    # Ensure start and end coordinates exist; if from task2, calculate from search_start
    first_exon = next(e for e in model['exons'] if e['is_first'])
    exon1_start = first_exon['start']
    exon1_end = first_exon['end']

    if 'guide_start' not in gene_guides.columns and 'protospacer_start0' in gene_guides.columns:
        search_start = max(1, exon1_start - 25)
        gene_guides['guide_start'] = gene_guides.apply(
            lambda r: search_start + int(r['protospacer_start0']) - (3 if str(r.get(strand_col, '+')).lower() in ('-', 'reverse', 'rev') else 0),
            axis=1
        )
        gene_guides['guide_end'] = gene_guides['guide_start'] + 22

    # Sort guides from lowest genomic start coordinate to highest
    gene_guides = gene_guides.sort_values('guide_start', ascending=True).reset_index(drop=True)
    num_guides = len(gene_guides)

    # Calculate target window bounds encompassing Exon 1 and all guides with clean padding
    min_coord = min([exon1_start] + gene_guides['guide_start'].tolist()) - 35
    max_coord = max([exon1_end] + gene_guides['guide_end'].tolist()) + 35

    # Dynamically scale figure height to provide generous vertical spacing for any pool size
    dyn_height = max(figsize[1], 4.5 + num_guides * 0.95)
    fig, (ax_macro, ax_micro) = plt.subplots(
        2, 1, figsize=(figsize[0], dyn_height), gridspec_kw={'height_ratios': [1, max(2.3, num_guides * 0.35)]}
    )
    plt.subplots_adjust(hspace=0.4)

    # -------------------------------------------------------------
    # PANEL A: Macro View - Full Gene Sequence Architecture
    # -------------------------------------------------------------
    ax_macro.set_title(
        f"A. Full Gene Sequence Architecture: {gene} on Chromosome {model['chrom']} ({model['chrom_accession']})\n"
        f"Locus: {model['gene_start']:,} - {model['gene_end']:,} bp | Strand: {model['strand']} | Primary Transcript: {model['transcript']} | System: {model['mode']}",
        fontsize=12, fontweight='bold', loc='left', pad=12
    )

    gene_span = model['gene_end'] - model['gene_start']
    ax_macro.plot(
        [model['gene_start'], model['gene_end']], [0, 0],
        color='#4A5568', lw=2.5, zorder=1
    )

    # Direction arrows indicating transcription strand orientation
    num_arrows = 9
    arrow_xs = np.linspace(model['gene_start'] + 0.08 * gene_span, model['gene_end'] - 0.08 * gene_span, num_arrows)
    arrow_symbol = '<' if model['strand'] == '-' else '>'
    for ax_pos in arrow_xs:
        ax_macro.text(
            ax_pos, 0, arrow_symbol, ha='center', va='center',
            fontsize=12, fontweight='bold', color='#2D3748', zorder=2
        )

    # Render exons along the gene body
    for exon in model['exons']:
        e_start = exon['start']
        e_end = exon['end']
        e_len = e_end - e_start + 1
        is_first = exon['is_first']

        color = '#E53E3E' if is_first else '#3182CE'
        height = 0.55 if is_first else 0.38
        y_pos = -height / 2

        rect = patches.FancyBboxPatch(
            (e_start, y_pos), max(e_len, 120), height,
            boxstyle="square,pad=0",
            facecolor=color, edgecolor='#1A202C', lw=1.2, zorder=3
        )
        ax_macro.add_patch(rect)

        # Highlight First Exon with dedicated pointer annotation
        if is_first:
            ax_macro.annotate(
                f"First Exon (Targeted)\n{e_start:,} - {e_end:,} ({e_len} bp)",
                xy=((e_start + e_end) / 2, height / 2),
                xytext=((e_start + e_end) / 2, 0.65),
                ha='center', va='bottom', fontsize=9.5, fontweight='bold',
                color='#C53030',
                arrowprops=dict(arrowstyle="->", color='#C53030', lw=1.8),
                zorder=5
            )
        else:
            if gene == 'CXCL11' or exon['exon_num'] in (2, 5, 9):
                ax_macro.text(
                    (e_start + e_end) / 2, -0.42, f"E{exon['exon_num']}",
                    ha='center', va='top', fontsize=8, color='#4A5568', fontweight='bold'
                )

    # Shaded CRISPR design target window
    zoom_box = patches.Rectangle(
        (min_coord, -0.65), max_coord - min_coord, 1.3,
        linewidth=1.8, edgecolor='#DD6B20', facecolor='#FEEBC8',
        alpha=0.45, linestyle='--', zorder=0
    )
    ax_macro.add_patch(zoom_box)

    ax_macro.set_xlim(model['gene_start'] - 0.03 * gene_span, model['gene_end'] + 0.03 * gene_span)
    ax_macro.set_ylim(-0.85, 1.1)
    ax_macro.set_yticks([])
    ax_macro.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
    ax_macro.set_xlabel("Genomic Coordinate (GRCh38, bp)", fontsize=10.5, fontweight='bold')
    ax_macro.grid(axis='x', linestyle=':', alpha=0.5)

    macro_legend_elements = [
        patches.Patch(facecolor='#E53E3E', edgecolor='#1A202C', label='First Exon (Targeted)'),
        patches.Patch(facecolor='#3182CE', edgecolor='#1A202C', label='Other Exons'),
        Line2D([0], [0], color='#4A5568', lw=2, label=f"Transcription ({model['strand']} strand)"),
        patches.Patch(facecolor='#FEEBC8', edgecolor='#DD6B20', linestyle='--', label='Target Window (Expanded Below)')
    ]
    ax_macro.legend(handles=macro_legend_elements, loc='upper right', frameon=True, fontsize=9, framealpha=0.95)

    # -------------------------------------------------------------
    # PANEL B: Micro View - Target Region, Exon 1 & gRNA / PAM Tracks
    # -------------------------------------------------------------
    ax_micro.set_title(
        f"B. Target Locus Architecture: First Exon & Selected gRNA / PAM Binding Sites\n"
        f"Exon 1 Window ({exon1_start:,} - {exon1_end:,}) with {num_guides} eligible CRISPR guides",        fontsize=12, fontweight='bold', loc='left', pad=12
    )

    num_guides = len(gene_guides)
    y_exon = num_guides + 1.1
    y_tracks = list(range(num_guides, 0, -1))

    # Highlight Exon 1 vertical span across all guide tracks
    ax_micro.axvspan(
        exon1_start, exon1_end,
        facecolor='#FFF5F5', alpha=0.6, zorder=0
    )

    # Draw First Exon track at the top
    exon_box = patches.FancyBboxPatch(
        (exon1_start, y_exon - 0.22), exon1_end - exon1_start, 0.44,
        boxstyle="square,pad=0",
        facecolor='#FEB2B2', edgecolor='#E53E3E', lw=2, zorder=3
    )
    ax_micro.add_patch(exon_box)
    ax_micro.text(
        (exon1_start + exon1_end) / 2, y_exon,
        f"First Exon: {exon1_start:,} - {exon1_end:,} ({exon1_end - exon1_start + 1} bp)",
        ha='center', va='center', fontsize=10, fontweight='bold', color='#9B2C2C', zorder=4
    )

    # Flanking DNA lines
    ax_micro.plot([min_coord, exon1_start], [y_exon, y_exon], color='#A0AEC0', lw=2, linestyle='--', zorder=2)
    ax_micro.plot([exon1_end, max_coord], [y_exon, y_exon], color='#A0AEC0', lw=2, linestyle='--', zorder=2)

    # Boundary annotations at top of Exon 1
    ax_micro.text(exon1_start, y_exon + 0.35, f"Start: {exon1_start:,}",
                  ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#C53030')
    ax_micro.text(exon1_end, y_exon + 0.35, f"End: {exon1_end:,}",
                  ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#C53030')

    # Exon boundary lines across all tracks
    ax_micro.axvline(exon1_start, color='#E53E3E', linestyle='--', lw=1.2, alpha=0.75, zorder=1)
    ax_micro.axvline(exon1_end, color='#E53E3E', linestyle='--', lw=1.2, alpha=0.75, zorder=1)

    y_tick_labels = []
    y_tick_positions = []

    for idx, (_, row) in enumerate(gene_guides.iterrows()):
        y = y_tracks[idx]
        guide_id = row[id_col]
        target_seq = str(row['targetSeq']).strip().upper()
        g_start = int(row['guide_start'])
        g_end = int(row['guide_end'])
        strand = row[strand_col]
        strand_sign = '+' if strand in ('+', 'forward', 'FORW') else '-'

        if strand_sign == '+':
            proto_start, proto_end = g_start, g_start + 19
            pam_start, pam_end = g_start + 20, g_end
            cut_pos = g_start + 16.5
            proto_seq = target_seq[:20]
            pam_seq = target_seq[20:]
            arrow_sym = '→'
        else:
            pam_start, pam_end = g_start, g_start + 2
            proto_start, proto_end = g_start + 3, g_end
            cut_pos = g_start + 5.5
            proto_seq = target_seq[:20]
            pam_seq = target_seq[20:]
            arrow_sym = '←'

        # Baseline horizontal reference
        ax_micro.plot([min_coord, max_coord], [y, y], color='#EDF2F7', lw=1, zorder=1)

        # 1. Protospacer bar (20 nt)
        proto_width = proto_end - proto_start + 1
        proto_box = patches.Rectangle(
            (proto_start, y - 0.22), proto_width, 0.44,
            facecolor='#2B6CB0', edgecolor='#1A365D', lw=1.2, zorder=3
        )
        ax_micro.add_patch(proto_box)

        # 2. PAM region (3 nt) in contrasting bold amber
        pam_width = pam_end - pam_start + 1
        pam_box = patches.Rectangle(
            (pam_start, y - 0.22), pam_width, 0.44,
            facecolor='#DD6B20', edgecolor='#7B341E', lw=1.2, zorder=4
        )
        ax_micro.add_patch(pam_box)

        # 3. Arrow on the protospacer bar
        arrow_x = (proto_start + proto_end) / 2
        ax_micro.text(
            arrow_x, y, f"{arrow_sym} 20nt",
            ha='center', va='center', fontsize=8.5, fontweight='bold', color='white', zorder=5
        )

        # 4. Cut Site Marker (triangle below the bar to avoid any text collision)
        ax_micro.plot(
            [cut_pos, cut_pos], [y - 0.22, y + 0.22],
            color='#E53E3E', linestyle='-', lw=2.2, zorder=6
        )
        ax_micro.scatter(
            [cut_pos], [y - 0.28],
            marker='^', color='#E53E3E', s=55, zorder=7
        )

        # 5. Clean Sequence label above the bar with safe clearance
        seq_text = f"5'-{proto_seq} [{pam_seq}]-3'" if strand_sign == '+' else f"3'-[{pam_seq}] {proto_seq}-5'"
        ax_micro.text(
            (g_start + g_end) / 2, y + 0.35, seq_text,
            ha='center', va='bottom', fontsize=8.2, fontfamily='monospace',
            color='#1A202C', fontweight='bold', zorder=5
        )

        # 6. Structured metadata in Y-tick label (guarantees ZERO in-plot overlap)
        score_val = row.get('task2_composite', row.get('Composite_Score', None))
        cfd_val = row.get('cfdSpecScore', None)
        eff_val = row.get('Efficiency', None)

        info_lines = [
            f"{guide_id} ({strand_sign})"
        ]
        if pd.notna(score_val):
            info_lines.append(f"Composite: {float(score_val):.1f}")
        elif pd.notna(cfd_val) and pd.notna(eff_val):
            info_lines.append(f"CFD:{float(cfd_val):.0f} Eff:{float(eff_val):.0f}")

        y_tick_labels.append("\n".join(info_lines))
        y_tick_positions.append(y)

    # Y-axis setup
    y_tick_positions.append(y_exon)
    y_tick_labels.append("GENE MODEL\n(First Exon)")
    ax_micro.set_yticks(y_tick_positions)
    ax_micro.set_yticklabels(y_tick_labels, fontsize=8.5, fontweight='normal')
    ax_micro.set_ylim(0.2, y_exon + 0.9)

    ax_micro.set_xlim(min_coord, max_coord)
    ax_micro.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
    ax_micro.set_xlabel("Genomic Coordinate (GRCh38, bp)", fontsize=10.5, fontweight='bold')
    ax_micro.grid(axis='x', linestyle=':', alpha=0.5)

    # Micro Legend placed horizontally below the plot (guaranteed ZERO plot overlap)
    micro_legend_elements = [
        patches.Patch(facecolor='#FEB2B2', edgecolor='#E53E3E', label='First Exon (Targeted Region)'),
        patches.Patch(facecolor='#2B6CB0', edgecolor='#1A365D', label='gRNA Protospacer (20 nt)'),
        patches.Patch(facecolor='#DD6B20', edgecolor='#7B341E', label='PAM Region (3 nt, NGG)'),
        Line2D([0], [0], color='#E53E3E', lw=2, marker='^', markersize=7, label='SpCas9 Cleavage Site (Cut)'),
        Line2D([0], [0], color='#A0AEC0', linestyle='--', lw=1.8, label='Flanking DNA')
    ]
    ax_micro.legend(
        handles=micro_legend_elements,
        loc='upper center', bbox_to_anchor=(0.5, -0.16),
        ncol=5, frameon=True, fontsize=9, framealpha=0.95
    )

    plt.tight_layout()

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_p, dpi=200, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_all_chromosome_targets(guides_df=None, output_dir=None, show=True, prefix="task3"):
    """
    Convenience wrapper to plot both chromosome targets (Chr 4 / CXCL11 and Chr 18 / SERPINB2).

    Parameters
    ----------
    guides_df : pd.DataFrame, optional
        DataFrame containing guide metadata. If None, loaded from manifest.
    output_dir : str or Path, optional
        Directory where generated figures will be stored.
    show : bool, default True
        Whether to display figures.
    prefix : str, default 'task3'
        Prefix for output filenames.

    Returns
    -------
    dict
        Mapping of gene name to generated matplotlib Figure.
    """
    figs = {}
    for chrom, gene in [('4', 'CXCL11'), ('18', 'SERPINB2')]:
        out_p = None
        if output_dir:
            out_p = Path(output_dir) / f"{prefix}_{gene.lower()}_grna_architecture.png"
        figs[gene] = plot_gene_grna_architecture(
            gene_or_chrom=gene,
            guides_df=guides_df,
            output_path=out_p,
            show=show
        )
    return figs
