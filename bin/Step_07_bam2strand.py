#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
Infer strand specificity from a BAM file, generate a TSV summary,
and produce a pie chart.

Also computes Splices per Kilobase (splice/kb) as a complementary
genomic-DNA (gDNA) contamination QC proxy. In an in-house EV RNA-seq pilot
(EXODUS-M + miRNeasy Advanced, 400 uL plasma; N=5), the single no-DNase
control had splice/kb = 0.13, whereas the four TURBO DNase pilot libraries
had splice/kb values of 1.67, 0.63, 1.64 and 0.94. These empirical values
are exploratory QC references for residual genomic-DNA contribution, not an
absolute gDNA assay, optimized DNase protocol or universal cutoff.

Metric definition:
  splice/kb = (total splice junction crossings) / (n_unique_reads x avg_mapped_len / 1000)

Both the STAR log path and the BAM fallback path count splice *junction crossings*
(not spliced reads). A read with CIGAR 30M500N30M200N30M crosses 2 junctions
and contributes 2, consistent with STAR "Number of splices: Total".

When gDNA-corrected outputs are requested, EVscope can apply opposite-strand
subtraction as a heuristic sensitivity analysis. Splice/kb is reported as an
independent, complementary QC proxy and should be interpreted with other QC
summaries.

Reference: EVscope (https://www.biorxiv.org/content/10.1101/2025.06.24.660984v1)

The default test read number is 100,000,000.
"""

import matplotlib as mpl
mpl.use("pdf")
import subprocess
import sys
import os
import re
import argparse
import matplotlib.pyplot as plt

# Set font properties for publication quality
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

# Compiled regex for splice junction detection in CIGAR
_SPLICE_RE = re.compile(r'(\d+)N')
# Compiled regex for aligned match/mismatch operations in CIGAR
_MATCH_RE  = re.compile(r'(\d+)[M=X]')

# Exploratory in-house DNase pilot references for the splice/kb panel.
# These values are displayed as QC references only. They are not universal
# thresholds and do not constitute direct genomic DNA quantification.
NO_DNASE_SPLICE_PER_KB = 0.13
DNASE_TREATED_SPLICE_PER_KB_VALUES = [1.67, 0.63, 1.64, 0.94]
MEAN_DNASE_TREATED_SPLICE_PER_KB = sum(DNASE_TREATED_SPLICE_PER_KB_VALUES) / len(DNASE_TREATED_SPLICE_PER_KB_VALUES)


def _safe_float(value):
    """Return float(value), or None for NA/non-finite plotting inputs."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if x != x or x in (float('inf'), float('-inf')):
        return None
    return x


def clean_sample_name_from_bam_basename(base_name):
    """Convert EVscope BAM basename to display sample name for plotting."""
    suffixes = [
        '_STAR_umi_dedup_Aligned.sortedByCoord.out',
        '_Aligned.sortedByCoord.out',
        '_Aligned.sortedByCoord_umi_dedup.out',
    ]
    sample_name = base_name
    for suffix in suffixes:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break
    return sample_name

def plot_strandness_and_splice_qc(
    output_dir, base_name, lib_type, frac_fwd, frac_rev, frac_failed, spkb
):
    """Plot strandness pie chart plus a right-side splice/kb QC barplot."""
    f_failed = float(frac_failed)
    f_forward = float(frac_fwd)
    f_reverse = float(frac_rev)

    labels = ['Forward strand', 'Reverse strand', 'Undetermined strand']
    sizes  = [f_forward, f_reverse, f_failed]
    colors = ['#FF4500', '#FFA500', '#87CEEB']

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(8.27, 3.0), dpi=300,
        gridspec_kw={'width_ratios': [1.0, 1.15]}
    )

    # --- Left: original strandness pie chart ---
    wedges, _ = ax1.pie(
        sizes, startangle=90, colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 0.8},
        radius=1.08, center=(0, 0)
    )
    ax1.axis('equal')
    ax1.set_title('')

    legend_labels = [f"{lbl} ({sz*100:.1f}%)" for lbl, sz in zip(labels, sizes)]
    fig.legend(
        wedges, legend_labels, loc='lower center', bbox_to_anchor=(0.255, 0.090),
        frameon=False, prop={'size': 6.0}, ncol=1, labelspacing=0.28,
        handlelength=1.0, handletextpad=0.38, borderaxespad=0.0
    )

    # --- Right: splice/kb QC reference barplot ---
    sample_spkb = _safe_float(spkb)
    display_sample_name = clean_sample_name_from_bam_basename(base_name)
    bar_labels = [display_sample_name, 'Mean DNase pilot (N=4)', 'No-DNase pilot (N=1)']
    bar_values = [sample_spkb if sample_spkb is not None else 0.0,
                  MEAN_DNASE_TREATED_SPLICE_PER_KB,
                  NO_DNASE_SPLICE_PER_KB]
    bar_colors = ['#2171B5', '#1A9850', '#D73027']

    y_pos = [0.00, 0.40, 0.80]
    ax2.barh(
        y_pos, bar_values, color=bar_colors, edgecolor='white',
        linewidth=0.35, height=0.24
    )
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(bar_labels, fontsize=6.6)
    ax2.tick_params(axis='y', pad=3, length=0)
    ax2.set_xlabel('Splice/kb', fontsize=7.2, labelpad=5)
    ax2.set_title('')

    xmax = max(1.5, max(bar_values) * 1.22)
    dx = xmax * 0.018
    ax2.set_xlim(0, xmax)
    ax2.set_ylim(1.05, -0.25)  # sample at top
    ax2.grid(axis='x', linestyle=':', linewidth=0.4, alpha=0.35)
    ax2.set_axisbelow(True)

    for y, value, label, color in zip(y_pos, bar_values, bar_labels, bar_colors):
        if label == display_sample_name and sample_spkb is None:
            ax2.text(dx, y, 'NA', va='center', ha='left', fontsize=6.2,
                     color='#666666', fontweight='bold')
        else:
            ax2.text(value + dx, y, f'{value:.2f}', va='center', ha='left',
                     fontsize=6.2, color=color, fontweight='bold')

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)

    # Manual compact layout matching the Step 18 two-panel style.
    ax1.set_position([0.075, 0.260, 0.355, 0.610])
    ax2.set_position([0.610, 0.285, 0.350, 0.550])
    fig.text(0.252, 0.930, f'Strand Specificity\n({lib_type})',
             ha='center', va='top', fontsize=8.2)
    fig.text(0.785, 0.930, 'splice/kb gDNA QC proxy',
             ha='center', va='top', fontsize=8.2)

    pie_pdf = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.pdf")
    pie_png = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.png")
    fig.savefig(pie_pdf, format='pdf')
    fig.savefig(pie_png, format='png', dpi=300)
    plt.close(fig)
    print(f"Combined strandness + splice/kb QC chart saved as: {pie_pdf}, {pie_png}")


def count_cigar_junctions(cigar):
    """Count the number of N operations (splice junctions) in a CIGAR string."""
    return len(_SPLICE_RE.findall(cigar))


def cigar_mapped_len(cigar):
    """Sum of M (match/mismatch) operations in CIGAR string.

    This matches STAR 'Average mapped length': aligned M/=/X bases are counted.
    Soft-clipped (S) and intronic (N) bases are excluded, consistent with STAR
    reporting splice junction statistics from uniquely mapped reads only.

    Example: '19M170458N61M' -> 80 nt  (intron N excluded)
             '1S65M'         -> 65 nt  (soft-clip S excluded)
    """
    return sum(int(n) for n in _MATCH_RE.findall(cigar))


def compute_splice_per_kb(bam_file, star_log=None, sample_n=500000):
    """
    Compute splices per kilobase (splice/kb).

    Metric: splice/kb = n_splice_junctions / (n_unique_reads x avg_mapped_len / 1000)

    Both paths count splice *junction crossings* (consistent with STAR
    "Number of splices: Total"). A read spanning 2 introns = 2 junctions.

    Primary path: STAR Log.final.out (fast, exact). For paired-end libraries,
    STAR reports uniquely mapped reads at the fragment/read-pair level, and
    Average mapped length corresponds to the total aligned bases from both
    mates per uniquely mapped fragment (not insert size).
    Fallback: samtools view with up to sample_n uniquely mapped alignment
    records (-q 255). The fallback is an approximate alignment-record estimate.

    Returns: (splice_per_kb, n_junctions, n_reads, avg_len, source) or None on failure.
    """
    if star_log and os.path.isfile(star_log):
        try:
            metrics = {}
            with open(star_log) as f:
                for line in f:
                    if "|" in line:
                        k, v = [x.strip() for x in line.split("|")]
                        metrics[k] = v.replace('%', '').strip()
            n_junctions = int(metrics.get("Number of splices: Total", 0))
            n_unique    = int(metrics.get("Uniquely mapped reads number", 0))
            avg_len     = float(metrics.get("Average mapped length", 0))
            if n_unique > 0 and avg_len > 0:
                splice_per_kb = n_junctions * 1000.0 / (n_unique * avg_len)
                return round(splice_per_kb, 4), n_junctions, n_unique, round(avg_len, 1), 'STAR_Log.final.out'
        except Exception as e:
            print(f"Warning: STAR log parse failed ({e}), falling back to BAM sampling.")

    # Fallback: sample uniquely mapped reads from BAM (-q 255 = unique in STAR)
    try:
        # Exclude secondary (0x100) and supplementary (0x800) alignments.
        # -q 255 is STAR-specific for uniquely mapped reads.
        cmd_total = [
            "samtools", "view", "-c", "-F", "2304", "-q", "255", bam_file
        ]
        r_total = subprocess.run(cmd_total, capture_output=True, text=True, check=True)
        n_total = int(r_total.stdout.strip())
        if n_total == 0:
            return 0.0, 0, 0, 0.0, 'BAM_alignment_record_estimate'

        frac = min(1.0, sample_n / n_total)
        sample_flag = ["-s", f"{frac:.6f}"] if frac < 1.0 else []
        # -q 255: uniquely mapped only for STAR output (consistent with STAR path)
        cmd_view = ["samtools", "view", "-F", "2304", "-q", "255"] + sample_flag + [bam_file]
        r_view = subprocess.run(cmd_view, capture_output=True, text=True, check=True)

        n_junctions_sampled = 0
        total_len = 0
        n_sampled = 0
        for line in r_view.stdout.splitlines():
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            n_sampled += 1
            n_junctions_sampled += count_cigar_junctions(fields[5])
            total_len += cigar_mapped_len(fields[5])  # aligned M/=/X ops, matches STAR avg_mapped_len

        if n_sampled == 0:
            return 0.0, 0, n_total, 0.0, 'BAM_alignment_record_estimate'

        avg_len = total_len / n_sampled
        junction_rate = n_junctions_sampled / n_sampled
        splice_per_kb = junction_rate * 1000.0 / avg_len if avg_len > 0 else 0.0
        n_junctions_est = int(junction_rate * n_total)
        return round(splice_per_kb, 4), n_junctions_est, n_total, round(avg_len, 1), 'BAM_alignment_record_estimate'

    except Exception as e:
        print(f"Warning: BAM-based splice/kb computation failed: {e}")
        return None


def run_infer_experiment(bam_file, refgene_bed, test_read_num, output_dir, star_log=None):
    """
    Runs RSeQC infer_experiment.py, parses the output, and generates a TSV
    and pie chart. Also computes and appends splice/kb to the TSV.
    """
    if not os.path.isfile(bam_file):
        print(f"Error: Input BAM file not found at {bam_file}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(bam_file).replace('.bam', '')
    tsv_path  = os.path.join(output_dir, f"{base_name}_bam2strandness.tsv")

    frac_fwd, frac_rev, frac_failed, lib_type = "0", "0", "0", "Unknown"
    spkb, n_junc, n_tot, avg_l, spkb_source = "NA", "NA", "NA", "NA", "NA"

    with open(tsv_path, "w") as f:
        f.write(
            "BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t"
            "%reads failed to determine:\t"
            "Splice_per_kb\tn_splice_junctions\tn_unique_reads\tavg_mapped_len_nt\tsplice_per_kb_source\n"
        )
        try:
            cmd  = ["infer_experiment.py", "-i", bam_file, "-r", refgene_bed,
                    "-s", str(test_read_num)]
            proc = subprocess.run(cmd, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            lines       = proc.stdout.strip().split('\n')
            lib_type    = lines[0].strip() if lines else "Library type not found"
            frac_failed = lines[1].split(": ")[1] if len(lines) > 1 else "0"
            frac_fwd    = lines[2].split(": ")[1] if len(lines) > 2 else "0"
            frac_rev    = lines[3].split(": ")[1] if len(lines) > 3 else "0"

            result = compute_splice_per_kb(bam_file, star_log=star_log)
            if result:
                spkb, n_junc, n_tot, avg_l, spkb_source = result

            f.write(
                f"{bam_file}\t{lib_type}\t{frac_fwd}\t{frac_rev}\t{frac_failed}\t"
                f"{spkb}\t{n_junc}\t{n_tot}\t{avg_l}\t{spkb_source}\n"
            )

        except subprocess.CalledProcessError as e:
            error_message = f"Error processing {bam_file}: {e.stderr}"
            f.write(f"{bam_file}\tError\tNA\tNA\tNA\tNA\tNA\tNA\tNA\terror\n")
            print(error_message)
            sys.exit(1)

    print(f"Strandness results written to {tsv_path}")
    print(f"Splice/kb = {spkb}  (n_junctions={n_junc}, n_reads={n_tot}, avg_len={avg_l}nt, source={spkb_source})")

    # --- Combined chart: strandness pie + right-side splice/kb QC barplot ---
    try:
        plot_strandness_and_splice_qc(
            output_dir, base_name, lib_type, frac_fwd, frac_rev, frac_failed, spkb
        )
    except ValueError:
        print("Error: Could not convert strand fractions to numbers. Skipping chart generation.")
        return


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Infer strand specificity and generate a TSV summary and pie chart. "
            "Also computes splices per kilobase (splice/kb) as a complementary "
            "gDNA-contamination QC proxy."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--input_bam",     required=True)
    parser.add_argument("--bed",           required=True,
                        help="Non-overlapping exon BED for infer_experiment.py.")
    parser.add_argument("--test_read_num", type=int, default=100000000)
    parser.add_argument("--output_dir",    required=True)
    parser.add_argument("--star_log",      required=False, default=None,
                        help="STAR Log.final.out from Step 6 refined alignment. "
                             "Primary source for splice/kb (fast, exact). "
                             "Falls back to BAM sampling if absent.")
    args = parser.parse_args()
    run_infer_experiment(
        args.input_bam, args.bed, args.test_read_num,
        args.output_dir, star_log=args.star_log
    )


if __name__ == "__main__":
    main()
