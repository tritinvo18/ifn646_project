import re
from pathlib import Path
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

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
        gene_guides['guide_start'] = search_start + gene_guides['protospacer_start0']
        gene_guides['guide_end'] = gene_guides['guide_start'] + 22

    # Sort guides from lowest genomic start coordinate to highest
    gene_guides = gene_guides.sort_values('guide_start', ascending=True).reset_index(drop=True)

    # Calculate target window bounds encompassing Exon 1 and all guides with clean padding
    min_coord = min([exon1_start] + gene_guides['guide_start'].tolist()) - 35
    max_coord = max([exon1_end] + gene_guides['guide_end'].tolist()) + 35

    # Create figure with 2 panels: Macro (full gene) and Micro (Exon 1 + gRNAs)
    fig, (ax_macro, ax_micro) = plt.subplots(
        2, 1, figsize=figsize, gridspec_kw={'height_ratios': [1, 2.3]}
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
        f"Exon 1 Window ({exon1_start:,} - {exon1_end:,}) with 4 Top-Ranked CRISPR Guides (Zero Overlap)",
        fontsize=12, fontweight='bold', loc='left', pad=12
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
            f"Guide: {guide_id} ({strand_sign})",
            f"PAM: {pam_seq} ({'3\'' if strand_sign == '+' else '5\''})",
            f"Coords: {g_start:,}..{g_end:,}"
        ]
        if pd.notna(score_val):
            info_lines.append(f"Score: {float(score_val):.1f}")
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


def plot_all_chromosome_targets(guides_df=None, output_dir=None, show=True):
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

    Returns
    -------
    dict
        Mapping of gene name to generated matplotlib Figure.
    """
    figs = {}
    for chrom, gene in [('4', 'CXCL11'), ('18', 'SERPINB2')]:
        out_p = None
        if output_dir:
            out_p = Path(output_dir) / f"task3_{gene.lower()}_grna_architecture.png"
        figs[gene] = plot_gene_grna_architecture(
            gene_or_chrom=gene,
            guides_df=guides_df,
            output_path=out_p,
            show=show
        )
    return figs
