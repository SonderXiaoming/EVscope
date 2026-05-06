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
import matplotlib.pyplot as plt
import numpy as np

# Publication-quality parameters
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


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

# Enrichment bar chart: sorted by enrichment (high to low)
# Colors: red/orange for high enrichment (Exon, 5'UTR), light blue for others
ENRICH_COLORS = {
    'Exon': '#E41A1C',
    "5'UTR": '#FF7F00',
    "3'UTR": '#6BAED6',
    'Promoter': '#6BAED6',
    'Downstream 2kb': '#6BAED6',
    'Intron': '#6BAED6',
    'Intergenic': '#6BAED6',
}


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

    # ── Combined figure: pie (left) + enrichment bar (right) ──
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5), dpi=300,
                                    gridspec_kw={'width_ratios': [1, 1.2]})

    # --- Left: Pie chart ---
    pie_counts = [region_counts[r] for r in DISPLAY_ORDER]
    pie_sizes = [c / total_reads for c in pie_counts]
    pie_colors = ['red', 'orange', 'cyan', 'green', 'purple', 'blue', 'pink']

    wedges, _ = ax1.pie(
        pie_sizes, startangle=90, colors=pie_colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        radius=1.0, center=(0, 0))
    ax1.axis('equal')
    ax1.set_title(f"Read Distribution\n({args.sampleName})", fontsize=8, pad=8)

    legend_labels = [f"{lbl} ({s*100:.1f}%)" for lbl, s in zip(DISPLAY_ORDER, pie_sizes)]
    ax1.legend(wedges, legend_labels, loc='center left',
               bbox_to_anchor=(0.85, 0.5), frameon=False, prop={'size': 5.5})

    # --- Right: Enrichment bar chart ---
    # Sort by enrichment score (high to low)
    sorted_regions = sorted(enrichment.keys(), key=lambda r: enrichment[r], reverse=True)
    sorted_scores = [enrichment[r] for r in sorted_regions]
    sorted_colors = [ENRICH_COLORS.get(r, '#6BAED6') for r in sorted_regions]

    bars = ax2.barh(range(len(sorted_regions)), sorted_scores,
                     color=sorted_colors, edgecolor='white', height=0.6)
    ax2.set_yticks(range(len(sorted_regions)))
    ax2.set_yticklabels(sorted_regions, fontsize=6)
    ax2.set_xlabel('Enrichment Score\n(observed / expected)', fontsize=7)
    ax2.set_title('Genomic Region Enrichment', fontsize=8)
    ax2.invert_yaxis()  # Highest at top

    # Annotate enrichment values
    for i, (score, region) in enumerate(zip(sorted_scores, sorted_regions)):
        ax2.text(score + 0.3, i, f'{score:.1f}', va='center', ha='left',
                 fontsize=6, color='#2171B5', fontweight='bold')

    # Add dashed line at enrichment = 1.0 (expected)
    ax2.axvline(x=1.0, color='grey', linestyle='--', linewidth=0.5, alpha=0.7)
    ax2.text(1.05, len(sorted_regions) - 0.3, 'expected', fontsize=5,
             color='grey', alpha=0.7)

    # Clean up spines
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.pdf")
    png_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.png")
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', pad_inches=0.02)
    fig.savefig(png_path, format='png', bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

    print(f"Plot saved: {png_path} and {pdf_path}")
