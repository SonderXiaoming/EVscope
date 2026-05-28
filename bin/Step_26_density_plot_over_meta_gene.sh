#!/bin/bash
# SPDX-License-Identifier: MIT
# Script to generate stacked profile plots over meta-gene regions from a bigWig file and BED files.
# The Y-axis scale is fixed (--yMin 0 --yMax 5) for consistent scaling across samples.
# Outputs PNG and SVG plots as ${output_prefix}_bed_stacked_profile_meta_gene.png/.svg.

# Function to display usage
usage() {
    echo "Usage: $0 --input_bw_file <bigWig_file> --input_bed_files \"[<bed_file1>,<bed_file2>,...]\" --input_bed_labels \"[<label1>,<label2>,...]\" --output_dir <output_directory> [--random_tested_row_num_per_bed <all|number>] [--blackListFileName <blacklist_bed_file>]"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_bw_file) bw_file="$2"; shift ;;
        --input_bed_files) input_bed_files="$2"; shift ;;
        --input_bed_labels) input_bed_labels="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        --random_tested_row_num_per_bed) random_tested_row_num_per_bed="$2"; shift ;;
        --threads) num_threads="$2"; shift ;;
        --blackListFileName) blackListFileName="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$bw_file" || -z "$input_bed_files" || -z "$input_bed_labels" || -z "$output_dir" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Default random_tested_row_num_per_bed to 'all' if not provided
random_tested_row_num_per_bed=${random_tested_row_num_per_bed:-all}

PYTHON_BIN="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
if [[ -z "$PYTHON_BIN" ]]; then
    echo "Error: Step 26 requires python3 or python in PATH."
    exit 1
fi
for required_tool in computeMatrix plotProfile; do
    if ! command -v "$required_tool" >/dev/null 2>&1; then
        echo "Error: Step 26 requires $required_tool in PATH."
        exit 1
    fi
done

# Convert JSON arrays to Bash arrays without splitting paths on spaces.
parse_json_array() {
    local json_input="$1"
    "$PYTHON_BIN" - "$json_input" <<'PYJSON'
import json
import sys
try:
    values = json.loads(sys.argv[1])
except Exception as exc:
    raise SystemExit(f"Invalid JSON array: {exc}")
if not isinstance(values, list):
    raise SystemExit("Expected a JSON array")
for value in values:
    print(str(value))
PYJSON
}
mapfile -t bed_files < <(parse_json_array "$input_bed_files")
mapfile -t labels < <(parse_json_array "$input_bed_labels")

# Validate that bed_files and labels arrays have the same length
if [[ ${#bed_files[@]} -ne ${#labels[@]} ]]; then
    echo "Error: Number of BED files (${#bed_files[@]}) does not match number of labels (${#labels[@]})."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Define output prefix based on bigWig file name
output_prefix=$(basename "$bw_file" .bw)

# Process BED files for random sampling and validate non-empty files
processed_bed_files=()
valid_labels=()
for i in "${!bed_files[@]}"; do
    bed_file="${bed_files[i]}"
    output_bed="$output_dir/$(basename "${bed_file%.*}")_processed.bed"

    if [[ "$random_tested_row_num_per_bed" != "all" ]]; then
        # Get total number of lines in the BED file
        total_lines=$(wc -l < "$bed_file")
        # Use min(total_lines, random_tested_row_num_per_bed) for sampling
        sample_lines=$(( total_lines < random_tested_row_num_per_bed ? total_lines : random_tested_row_num_per_bed ))

        if [[ $sample_lines -gt 0 ]]; then
            # Randomly sample lines using shuf
            shuf "$bed_file" | head -n "$sample_lines" > "$output_bed"
        else
            # Skip empty BED files
            echo "Warning: BED file ${bed_file} is empty, skipping."
            continue
        fi
    else
        # Use the original BED file
        cp "$bed_file" "$output_bed"
    fi

    # Check if processed BED file is non-empty
    if [[ -s "$output_bed" ]]; then
        processed_bed_files+=("$output_bed")
        valid_labels+=("${labels[i]}")
    else
        echo "Warning: Processed BED file ${output_bed} is empty, skipping."
        rm -f "$output_bed"
    fi
done

# Check if there are valid BED files to process
if [[ ${#processed_bed_files[@]} -eq 0 ]]; then
    echo "Error: No valid non-empty BED files to process."
    exit 1
fi

# Set number of colors based on the number of valid BED files
num_colors=${#processed_bed_files[@]}

# Generate color list from matplotlib tab20 colormap for vibrant and distinct colors
python_code=$(cat <<EOF
import matplotlib.pyplot as plt
import numpy as np
import sys
num_colors = int(sys.argv[1])
# Use tab20 for up to 20 colors, otherwise interpolate using viridis for scalability
if num_colors <= 20:
    palette = plt.get_cmap("tab20").colors
    colors = palette[:num_colors]
else:
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0, 1, num_colors))
colors_hex = []
for color in colors:
    r, g, b = color[:3]
    colors_hex.append("#%02x%02x%02x" % (int(r*255), int(g*255), int(b*255)))
print(" ".join(colors_hex))
EOF
)
if ! color_list=$("$PYTHON_BIN" -c "$python_code" "$num_colors"); then
    echo "Error: Failed to generate Step 26 colors. Ensure matplotlib and numpy are installed for $PYTHON_BIN."
    exit 1
fi
IFS=' ' read -r -a colors <<< "$color_list"

# Check if computeMatrix output exists
matrix_file="$output_dir/${output_prefix}_bed_stacked_matrix_meta_gene.gz"
sorted_regions_file="$output_dir/${output_prefix}_bed_stacked_sorted_regions_meta_gene.bed"

if [ ! -f "$matrix_file" ]; then
    echo "computeMatrix output not found, running computeMatrix..."
    computeMatrix scale-regions \
        -S "$bw_file" \
        -R "${processed_bed_files[@]}" \
        --beforeRegionStartLength 0 \
        --regionBodyLength 100 \
        --afterRegionStartLength 0 \
        -o "$matrix_file" \
        --binSize 5 \
        -p ${num_threads:-4} \
        --missingDataAsZero \
        ${blackListFileName:+--blackListFileName "$blackListFileName"} \
        --outFileSortedRegions "$sorted_regions_file"
else
    echo "computeMatrix output exists, skipping computeMatrix..."
fi

# Create region labels with counts from sorted regions
regions_label=()
if [ -f "$sorted_regions_file" ]; then
    for i in "${!processed_bed_files[@]}"; do
        count=$(grep -c "$(basename "${processed_bed_files[i]}")" "$sorted_regions_file" 2>/dev/null || echo 0)
        regions_label+=("${valid_labels[i]} (${count})")
    done
else
    echo "Warning: Sorted regions file not found. Label counts will be set to zero."
    for i in "${!processed_bed_files[@]}"; do
        regions_label+=("${valid_labels[i]} (0)")
    done
fi

# Run plotProfile to generate PNG output
echo "Running plotProfile for PNG..."
plotProfile -m "$matrix_file" \
    -out "$output_dir/${output_prefix}_bed_stacked_profile_meta_gene.png" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 10 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Reads coverage" \
    --yMin 0 --yMax 5 \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}" || {
        echo "Error: plotProfile failed for PNG output."
        exit 1
    }

# Verify PNG file exists
if [ -f "$output_dir/${output_prefix}_bed_stacked_profile_meta_gene.png" ]; then
    echo "PNG profile generated at $output_dir/${output_prefix}_bed_stacked_profile_meta_gene.png"
else
    echo "Error: PNG file was not generated."
    exit 1
fi

# Run plotProfile to generate SVG output
echo "Running plotProfile for SVG..."
plotProfile -m "$matrix_file" \
    -out "$output_dir/${output_prefix}_bed_stacked_profile_meta_gene.svg" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 10 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Reads coverage" \
    --yMin 0 --yMax 5 \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}" || {
        echo "Error: plotProfile failed for SVG output."
        exit 1
    }

# Verify SVG file exists
if [ -f "$output_dir/${output_prefix}_bed_stacked_profile_meta_gene.svg" ]; then
    echo "SVG profile generated at $output_dir/${output_prefix}_bed_stacked_profile_meta_gene.svg"
else
    echo "Error: SVG file was not generated."
    exit 1
fi

# Clean up processed BED files
for bed_file in "${processed_bed_files[@]}"; do
    rm -f "$bed_file"
done
