# Example usage:
#   1) Set ANALYSIS_DIR to the analysis directory containing gene_matrix.txt
#      and condition_matrix.txt files.
#   2) Run this script to generate normalized matrices, heatmaps, and Excel
#      reports.

# Set paths to your files
ANALYSIS_DIR="/users/paxak6/Haploid_screening_ONT/barcode07_analysis_20250517_042959"
GENE_MATRIX="$ANALYSIS_DIR/condition_comparison/gene_matrix.txt"
CHROM_MATRIX="$ANALYSIS_DIR/condition_comparison/condition_matrix.txt"
OUTPUT_DIR="$ANALYSIS_DIR/matrix_visualization"

# Create output directory
mkdir -p $OUTPUT_DIR

# Run the script with both matrices using the enhanced version
# --top_genes: Number of genes to include in Excel files (200)
# --heatmap_genes: Number of genes to include in heatmaps (50)
python Normalize_Matrix_Visualization_Enhanced_v2.py \
    --matrix "$GENE_MATRIX" \
    --chrom_matrix "$CHROM_MATRIX" \
    --output_dir "$OUTPUT_DIR" \
    --top_genes 200 \
    --heatmap_genes 50 \
    --excel


