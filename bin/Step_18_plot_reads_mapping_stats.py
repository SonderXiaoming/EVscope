#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
Step 18: Compute read counts per genomic region from featureCounts TSV files.
Outputs:
  1. Summary TSV with read counts, fractions, and enrichment scores
  2. Combined figure: pie chart (left) + enrichment bar chart (right)
"""
import os
import sys
import argparse
import matplotlib as mpl
mpl.use("pdf")
import numpy as np

mpl.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 5,
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'lines.linewidth': 0.5,
})
import matplotlib.pyplot as plt


def sum_read_counts(tsv_file):
    """Sum counts in the last column of a TSV file, skipping headers and comments."""
    total = 0.0
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('Geneid'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    try:
                        total += float(parts[-1])
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        print(f"Warning: File not found {tsv_file}. Counting as 0.", file=sys.stderr)
        return 0.0
    return total


def compute_genome_fractions(saf_file):
    """Compute genome region sizes from merged SAF annotation file."""
    sizes = {}
    with open(saf_file, 'r') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            gene_id = parts[0]
            region = gene_id.split('__')[0]
            length = int(parts[3]) - int(parts[2])
            sizes[region] = sizes.get(region, 0) + length
    total = sum(sizes.values())
    fractions = {r: s / total for r, s in sizes.items()}
    return sizes, fractions


# Mapping from SAF region prefix to display label
SAF_TO_LABEL = {
    '5UTR': "5'UTR",
    'exon': 'Exon',
    '3UTR': "3'UTR",
    'intron': 'Intron',
    'promoter': 'Promoter',
    'downstream': 'Downstream 2kb',
    'intergenic': 'Intergenic',
}

# Display order for pie chart legend (same as original)
DISPLAY_ORDER = ["5'UTR", 'Exon', "3'UTR", 'Intron', 'Promoter', 'Downstream 2kb', 'Intergenic']

# Display order for the sample-level enrichment barplot.
BARPLOT_ORDER = ['Exon', "5'UTR", "3'UTR", 'Promoter', 'Downstream 2kb', 'Intron', 'Intergenic']

# Original read-distribution pie-chart palette. The enrichment panel reuses the
# same region colors so that each genomic feature has one visual identity.
REGION_COLORS = {
    "5'UTR": '#FF0000',        # original red
    'Exon': '#FFA500',           # original orange
    "3'UTR": '#00FFFF',        # original cyan
    'Intron': '#008000',         # original green
    'Promoter': '#800080',       # original purple
    'Downstream 2kb': '#000080', # original dark blue
    'Intergenic': '#F4B6C6',     # original light pink
}
ENRICH_COLORS = REGION_COLORS


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot read mapping statistics: pie chart + enrichment bar chart.'
    )
    parser.add_argument('--input_5UTR_readcounts', required=True)
    parser.add_argument('--input_exon_readcounts', required=True)
    parser.add_argument('--input_3UTR_readcounts', required=True)
    parser.add_argument('--input_intron_readcounts', required=True)
    parser.add_argument('--input_promoters_readcounts', required=True)
    parser.add_argument('--input_downstream_2Kb_readcounts', required=True)
    parser.add_argument('--input_intergenic_readcounts', required=True)
    parser.add_argument('--input_ENCODE_blacklist_readcounts', required=False, default=None,
                        help='(Deprecated) TSV for ENCODE blacklist counts')
    parser.add_argument('--sampleName', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--saf_file', required=False, default=None,
                        help='Merged SAF annotation file for enrichment calculation. '
                             'If not provided, uses built-in HG38 genome fractions.')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Sum reads per region ──
    region_counts = {
        "5'UTR":          sum_read_counts(args.input_5UTR_readcounts),
        'Exon':           sum_read_counts(args.input_exon_readcounts),
        "3'UTR":          sum_read_counts(args.input_3UTR_readcounts),
        'Intron':         sum_read_counts(args.input_intron_readcounts),
        'Promoter':       sum_read_counts(args.input_promoters_readcounts),
        'Downstream 2kb': sum_read_counts(args.input_downstream_2Kb_readcounts),
        'Intergenic':     sum_read_counts(args.input_intergenic_readcounts),
    }
    total_reads = sum(region_counts.values())
    if total_reads == 0:
        print('Error: Total read count is zero.', file=sys.stderr)
        sys.exit(1)

    read_fractions = {r: c / total_reads for r, c in region_counts.items()}

    # ── Compute genome fractions (for enrichment) ──
    if args.saf_file and os.path.exists(args.saf_file):
        genome_sizes, genome_fractions_raw = compute_genome_fractions(args.saf_file)
        # Map SAF keys to display labels
        genome_fractions = {}
        for saf_key, label in SAF_TO_LABEL.items():
            if saf_key in genome_fractions_raw:
                genome_fractions[label] = genome_fractions_raw[saf_key]
    else:
        # Built-in HG38 genome fractions (from HG38_all7_metagene_noOverlap.saf)
        genome_fractions = {
            "5'UTR":          0.003372,
            'Exon':           0.031191,
            "3'UTR":          0.016072,
            'Intron':         0.525170,
            'Promoter':       0.017601,
            'Downstream 2kb': 0.016987,
            'Intergenic':     0.389607,
        }

    # ── Enrichment scores ──
    enrichment = {}
    for region in DISPLAY_ORDER:
        obs = read_fractions.get(region, 0)
        exp = genome_fractions.get(region, 0.001)
        enrichment[region] = obs / exp

    # ── Save summary TSV ──
    tsv_path = os.path.join(args.output_dir,
                            f"{args.sampleName}_reads_mapping_readcounts_summary.tsv")
    with open(tsv_path, 'w') as f:
        f.write("Genomic_feature\tRead_counts\tRead_fraction(%)\tGenome_fraction(%)\tEnrichment_Score\n")
        for region in DISPLAY_ORDER:
            rc = int(region_counts.get(region, 0))
            rf = read_fractions.get(region, 0) * 100
            gf = genome_fractions.get(region, 0) * 100
            es = enrichment.get(region, 0)
            f.write(f"{region}\t{rc}\t{rf:.2f}\t{gf:.2f}\t{es:.2f}\n")
    print(f"Summary table saved: {tsv_path}")

    # ── Combined figure: original read-distribution pie (left) + enrichment panel (right) ──
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(6.6, 3.2), dpi=300,
        gridspec_kw={'width_ratios': [1.12, 1.0]}
    )

    # --- Left: Pie chart; keep original visual style ---
    pie_counts = [region_counts[r] for r in DISPLAY_ORDER]
    pie_sizes = [c / total_reads for c in pie_counts]
    pie_colors = [REGION_COLORS[r] for r in DISPLAY_ORDER]

    wedges, _ = ax1.pie(
        pie_sizes, startangle=90, colors=pie_colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        radius=1.05, center=(0, 0))
    ax1.axis('equal')
    ax1.set_title('')

    legend_labels = [f"{lbl} ({s*100:.1f}%)" for lbl, s in zip(DISPLAY_ORDER, pie_sizes)]
    fig.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(0.128, 0.455),
               frameon=False, prop={'size': 5.4}, ncol=1, labelspacing=0.32,
               handlelength=1.0, handletextpad=0.38, borderaxespad=0.0)

    # --- Right: Enrichment bar chart ---
    # Sort by enrichment score (high to low)
    sorted_regions = sorted(enrichment.keys(), key=lambda r: enrichment[r], reverse=True)
    sorted_scores = [enrichment[r] for r in sorted_regions]
    sorted_colors = [REGION_COLORS.get(r, '#999999') for r in sorted_regions]

    y_step = 0.33
    bar_height = 0.20
    y_pos = np.arange(len(sorted_regions)) * y_step
    y_mid = (y_pos[0] + y_pos[-1]) / 2
    xmax = max(1.5, max(sorted_scores) * 1.24)
    dx = xmax * 0.018

    bars = ax2.barh(y_pos, sorted_scores,
                     color=sorted_colors, edgecolor='white', linewidth=0.25, height=bar_height)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(sorted_regions, fontsize=5.5)
    ax2.tick_params(axis='y', pad=2)
    ax2.set_xlabel('Enrichment Score (fold)\nobserved / expected', fontsize=7, labelpad=4)
    ax2.set_title('')
    ax2.set_xlim(0, xmax)
    # Match the enrichment panel's visual height to the pie chart diameter.
    ax2.set_ylim(y_mid + 1.10, y_mid - 1.10)  # Highest enrichment at top

    # Annotate enrichment values
    for y, score, region in zip(y_pos, sorted_scores, sorted_regions):
        ax2.text(score + dx, y, f'{score:.1f}', va='center', ha='left',
                 fontsize=5.5, color='#2171B5', fontweight='bold')

    # Add dashed line at enrichment = 1.0 (expected)
    ax2.axvline(x=1.0, color='grey', linestyle='--', linewidth=0.5, alpha=0.7)
    ax2.text(1.0, -0.035, 'expected', transform=ax2.get_xaxis_transform(),
             rotation=-90, ha='center', va='top', fontsize=5.2,
             color='#666666', alpha=0.9)

    # Clean up spines
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)

    # Manual compact landscape layout: legend | pie | enrichment.
    # Titles are figure-level text so both panels start at exactly the same height.
    ax1.set_position([0.265, 0.205, 0.300, 0.620])
    ax2.set_position([0.645, 0.235, 0.320, 0.560])
    fig.text(0.415, 0.910, f'Read Distribution\n({args.sampleName})',
             ha='center', va='top', fontsize=8)
    fig.text(0.805, 0.910, f'Genomic Region Enrichment\n({args.sampleName})',
             ha='center', va='top', fontsize=8)

    # Save
    pdf_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.pdf")
    png_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.png")
    fig.savefig(pdf_path, format='pdf')
    fig.savefig(png_path, format='png', dpi=300)
    plt.close(fig)

    print(f"Plot saved: {png_path} and {pdf_path}")
