#!/usr/bin/env python3
"""
Normalize_Matrix_Visualization_Enhanced_v2.py
--------------------------------
Creates normalized and improved visualizations for condition comparison matrices
from haploid screening data. Normalizes counts by condition total to account for
differences in sequencing depth.

This script handles both:
1. Gene matrix (gene_matrix.txt) for gene-level analysis
2. Condition matrix (condition_matrix.txt) for chromosome-level analysis

Dependencies:
- pandas
- matplotlib
- seaborn
- numpy
- xlsxwriter (for Excel export)

Usage Example:
  python Normalize_Matrix_Visualization_Enhanced_v2.py \
    --matrix gene_matrix.txt \
    --output_dir matrix_viz \
    --top_genes 50 \
    --chrom_matrix condition_matrix.txt
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import datetime

# For Excel export
try:
    import xlsxwriter
    HAS_XLSXWRITER = True
except ImportError:
    HAS_XLSXWRITER = False

# Try to import seaborn for better visualizations
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}", file=sys.stderr)

def load_matrix(matrix_file):
    """
    Load a gene-condition matrix from a TSV file
    
    Args:
        matrix_file (str): Path to the matrix file
        
    Returns:
        pandas.DataFrame: The loaded matrix
    """
    try:
        # Read the matrix file
        df = pd.read_csv(matrix_file, sep='\t')
        
        # Check if required columns exist
        required_cols = ["Gene", "GeneID"]
        if not all(col in df.columns for col in required_cols):
            log(f"ERROR: Matrix file is missing required columns: {required_cols}")
            return None
        
        log(f"Loaded matrix with {df.shape[0]} genes and {df.shape[1]} columns")
        return df
    
    except Exception as e:
        log(f"ERROR loading matrix file: {str(e)}")
        return None

def load_chromosome_matrix(chrom_matrix_file):
    """
    Load a chromosome-condition matrix from a TSV file
    
    Args:
        chrom_matrix_file (str): Path to the chromosome matrix file
        
    Returns:
        pandas.DataFrame: The loaded chromosome matrix
    """
    try:
        # Read the matrix file with special handling for problematic formatting
        # Fix issues with missing tabs and newlines in the data
        with open(chrom_matrix_file, 'r') as f:
            lines = []
            for line in f:
                # Replace multiple whitespace with tabs if needed
                line = line.strip()
                if line:
                    lines.append(line)
        
        # Join lines and manually parse as TSV
        df = pd.read_csv(chrom_matrix_file, sep='\t')
        
        # Clean up the dataframe
        # Remove rows with missing values if needed
        df = df.dropna(how='all')
        
        # Clean columns - ensures all values are numeric
        for col in df.columns:
            if col != 'Chromosome':
                df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
                
        log(f"Loaded chromosome matrix with {df.shape[0]} chromosomes and {df.shape[1]-1} conditions")
        return df
    
    except Exception as e:
        log(f"ERROR loading chromosome matrix file: {str(e)}")
        log(f"Exception details: {e}")
        return None

def create_normalized_matrix(df):
    """
    Create a normalized version of the matrix
    
    Args:
        df (pandas.DataFrame): Original matrix
        
    Returns:
        pandas.DataFrame: Normalized matrix
    """
    log("Creating normalized matrix...")
    
    # Create a copy of the DataFrame
    norm_df = df.copy()
    
    # Get condition columns (excluding metadata columns)
    metadata_cols = ["Gene", "GeneID", "Total"]
    condition_cols = [col for col in df.columns if col not in metadata_cols]
    
    if not condition_cols:
        log("ERROR: No condition columns found in the matrix")
        return None
    
    # Calculate column totals for normalization
    column_totals = {}
    for col in condition_cols:
        column_totals[col] = df[col].sum()
        if column_totals[col] == 0:
            log(f"WARNING: Condition {col} has zero total counts")
            column_totals[col] = 1  # Avoid division by zero
    
    # Normalize each condition column
    for col in condition_cols:
        norm_df[f"{col}_normalized"] = (df[col] / column_totals[col]) * 100
    
    log(f"Created normalized matrix with {len(condition_cols)} additional columns")
    return norm_df

def create_normalized_chromosome_matrix(df):
    """
    Create a normalized version of the chromosome matrix
    
    Args:
        df (pandas.DataFrame): Original chromosome matrix
        
    Returns:
        pandas.DataFrame: Normalized chromosome matrix
    """
    log("Creating normalized chromosome matrix...")
    
    # Create a copy of the DataFrame
    norm_df = df.copy()
    
    # Get condition columns (excluding Chromosome column)
    condition_cols = [col for col in df.columns if col != 'Chromosome']
    
    if not condition_cols:
        log("ERROR: No condition columns found in the chromosome matrix")
        return None
    
    # Calculate column totals for normalization
    column_totals = {}
    for col in condition_cols:
        column_totals[col] = df[col].sum()
        if column_totals[col] == 0:
            log(f"WARNING: Condition {col} has zero total counts")
            column_totals[col] = 1  # Avoid division by zero
    
    # Normalize each condition column
    for col in condition_cols:
        norm_df[f"{col}_normalized"] = (df[col] / column_totals[col]) * 100
    
    log(f"Created normalized chromosome matrix with {len(condition_cols)} additional columns")
    return norm_df

def create_visualizations(df, norm_df, output_dir, top_n=50, heatmap_genes=50):
    """
    Create matrix visualizations
    
    Args:
        df (pandas.DataFrame): Original matrix
        norm_df (pandas.DataFrame): Normalized matrix
        output_dir (str): Output directory for visualizations
        top_n (int): Number of top genes to include in general visualizations
        heatmap_genes (int): Number of top genes to include in heatmaps
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get condition columns (excluding metadata columns)
    metadata_cols = ["Gene", "GeneID", "Total"]
    condition_cols = [col for col in df.columns if col not in metadata_cols]
    normalized_cols = [f"{col}_normalized" for col in condition_cols]
    
    if not condition_cols:
        log("ERROR: No condition columns found for visualization")
        return
    
    log(f"Creating visualizations for {len(condition_cols)} conditions...")
    
    # 1. Raw counts heatmap for top genes
    try:
        create_heatmap(
            df, condition_cols, 
            os.path.join(output_dir, "top_genes_heatmap_raw.png"),
            "Raw Insertion Counts Across Conditions",
            top_n=heatmap_genes  # Use heatmap_genes parameter instead of top_n
        )
    except Exception as e:
        log(f"ERROR creating raw counts heatmap: {str(e)}")
    
    # 2. Normalized counts heatmap for top genes
    try:
        create_heatmap(
            norm_df, normalized_cols, 
            os.path.join(output_dir, "top_genes_heatmap_normalized.png"),
            "Normalized Insertion Counts (% of condition total)",
            top_n=heatmap_genes  # Use heatmap_genes parameter instead of top_n
        )
    except Exception as e:
        log(f"ERROR creating normalized heatmap: {str(e)}")
    
    # 3. Condition totals bar chart
    try:
        create_condition_totals_chart(
            df, condition_cols,
            os.path.join(output_dir, "condition_totals.png")
        )
    except Exception as e:
        log(f"ERROR creating condition totals chart: {str(e)}")
    
    # 4. Distribution of insertions per gene
    try:
        create_insertion_distribution(
            df,
            os.path.join(output_dir, "insertion_distribution.png")
        )
    except Exception as e:
        log(f"ERROR creating insertion distribution: {str(e)}")
    
    # 5. Gene counts across conditions - limit to 10 genes for readability
    max_genes_for_comparison = min(10, heatmap_genes)
    try:
        create_gene_comparison_chart(
            df, condition_cols,
            os.path.join(output_dir, "top_genes_comparison.png"),
            top_n=max_genes_for_comparison
        )
    except Exception as e:
        log(f"ERROR creating gene comparison chart: {str(e)}")
    
    # 6. Normalized gene counts - limit to 10 genes for readability
    try:
        create_gene_comparison_chart(
            norm_df, normalized_cols,
            os.path.join(output_dir, "top_genes_comparison_normalized.png"),
            top_n=max_genes_for_comparison,
            normalized=True
        )
    except Exception as e:
        log(f"ERROR creating normalized gene comparison chart: {str(e)}")
    
    log(f"Created visualizations in {output_dir}")

def create_chromosome_visualizations(chrom_df, norm_chrom_df, output_dir):
    """
    Create chromosome-specific visualizations
    
    Args:
        chrom_df (pandas.DataFrame): Original chromosome matrix
        norm_chrom_df (pandas.DataFrame): Normalized chromosome matrix
        output_dir (str): Output directory for visualizations
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get condition columns
    condition_cols = [col for col in chrom_df.columns if col != 'Chromosome']
    normalized_cols = [f"{col}_normalized" for col in condition_cols]
    
    if not condition_cols:
        log("ERROR: No condition columns found for chromosome visualization")
        return
    
    log(f"Creating chromosome visualizations for {len(condition_cols)} conditions...")
    
    # 1. Filter to just standard chromosomes for clarity
    std_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                  '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                  '20', '21', '22', 'X', 'Y']
    
    # Create filtered datasets
    filtered_chrom_df = chrom_df[chrom_df['Chromosome'].isin(std_chroms)].copy()
    filtered_norm_df = norm_chrom_df[norm_chrom_df['Chromosome'].isin(std_chroms)].copy()
    
    # 2. Sort chromosomes naturally (1, 2, ... 10, 11)
    def chrom_sort_key(chrom_str):
        if chrom_str == 'X':
            return 23
        elif chrom_str == 'Y':
            return 24
        else:
            try:
                return int(chrom_str)
            except ValueError:
                return 999  # Place non-standard chroms at the end
    
    filtered_chrom_df['sort_key'] = filtered_chrom_df['Chromosome'].apply(chrom_sort_key)
    filtered_chrom_df = filtered_chrom_df.sort_values('sort_key')
    filtered_chrom_df = filtered_chrom_df.drop('sort_key', axis=1)
    
    filtered_norm_df['sort_key'] = filtered_norm_df['Chromosome'].apply(chrom_sort_key)
    filtered_norm_df = filtered_norm_df.sort_values('sort_key')
    filtered_norm_df = filtered_norm_df.drop('sort_key', axis=1)
    
    # 3. Chromosome insertion heatmap (raw)
    try:
        create_chromosome_heatmap(
            filtered_chrom_df, condition_cols,
            os.path.join(output_dir, "chromosome_heatmap_raw.png"),
            "Raw Insertion Counts by Chromosome",
        )
    except Exception as e:
        log(f"ERROR creating raw chromosome heatmap: {str(e)}")
    
    # 4. Chromosome insertion heatmap (normalized)
    try:
        create_chromosome_heatmap(
            filtered_norm_df, normalized_cols,
            os.path.join(output_dir, "chromosome_heatmap_normalized.png"),
            "Normalized Insertion Counts by Chromosome (% of condition total)",
        )
    except Exception as e:
        log(f"ERROR creating normalized chromosome heatmap: {str(e)}")
    
    # 5. Condition-specific insertion profiles
    try:
        create_chromosome_profile_chart(
            filtered_chrom_df, condition_cols,
            os.path.join(output_dir, "chromosome_profiles.png")
        )
    except Exception as e:
        log(f"ERROR creating chromosome profiles chart: {str(e)}")
    
    # 6. Normalized condition-specific profiles
    try:
        create_chromosome_profile_chart(
            filtered_norm_df, normalized_cols,
            os.path.join(output_dir, "chromosome_profiles_normalized.png"),
            normalized=True
        )
    except Exception as e:
        log(f"ERROR creating normalized chromosome profiles: {str(e)}")
    
    # 7. Treatment comparison of specific chromosome groups
    try:
        create_chromosome_group_comparison(
            filtered_chrom_df, filtered_norm_df, condition_cols,
            os.path.join(output_dir, "chromosome_group_comparison.png")
        )
    except Exception as e:
        log(f"ERROR creating chromosome group comparison: {str(e)}")
    
    log(f"Created chromosome visualizations in {output_dir}")

def save_figure_formats(plt, output_base, dpi=300):
    """
    Save figure in multiple formats (PNG and PDF)
    
    Args:
        plt: Matplotlib pyplot instance
        output_base (str): Base output path without extension
        dpi (int): DPI for the PNG output
    """
    # Save PNG version (standard)
    png_path = f"{output_base}.png"
    plt.savefig(png_path, dpi=dpi)
    
    # Save PDF version (A4 optimized)
    pdf_path = f"{output_base}.pdf"
    try:
        # A4 size in inches (8.27 × 11.69)
        # Use reduced target size to ensure margins (90% of A4)
        target_width = 8.27 * 0.85  # 85% of A4 width
        target_height = 11.69 * 0.85  # 85% of A4 height
        
        current_fig = plt.gcf()
        current_size = current_fig.get_size_inches()
        
        # Get the tight bbox to see how much margin space we need
        bbox_inches = plt.gcf().get_tightbbox(plt.gcf().canvas.get_renderer())
        
        # Add extra padding for margins (especially for labels and legends)
        # This ensures nothing is cut off at the edges
        plt.tight_layout(pad=1.5)
        
        # Keep aspect ratio while fitting to target width
        new_height = current_size[1] * (target_width / current_size[0])
        
        # If too tall for target height, scale down to fit
        if new_height > target_height:
            scale_factor = target_height / new_height
            new_height = target_height
            new_width = current_size[0] * scale_factor
        else:
            new_width = target_width
        
        # Temporarily resize for PDF save
        current_fig.set_size_inches(new_width, new_height)
        
        # Use bbox_inches='tight' to ensure all elements are included
        plt.savefig(pdf_path, format='pdf', dpi=dpi, bbox_inches='tight')
        
        # Restore original size
        current_fig.set_size_inches(current_size)
        
        return [png_path, pdf_path]
    except Exception as e:
        log(f"WARNING: Could not save PDF version: {str(e)}")
        return [png_path]

def create_heatmap(df, data_cols, output_file, title, top_n=50):
    """
    Create a heatmap visualization
    
    Args:
        df (pandas.DataFrame): Data frame
        data_cols (list): Data columns to include
        output_file (str): Output file path
        title (str): Chart title
        top_n (int): Number of top genes to include
    """
    # Sort by Total and take top N genes
    if "Total" in df.columns:
        top_df = df.sort_values("Total", ascending=False).head(top_n)
    else:
        # If no Total column, sum the data columns
        df["sum"] = df[data_cols].sum(axis=1)
        top_df = df.sort_values("sum", ascending=False).head(top_n)
    
    # Extract data for heatmap
    heatmap_data = top_df[data_cols].values
    genes = top_df["Gene"].values
    
    # Create figure with appropriate size
    plt.figure(figsize=(max(8, len(data_cols)), max(8, len(genes)/4)))
    
    # Use seaborn if available for better visualization
    if HAS_SEABORN:
        ax = sns.heatmap(
            heatmap_data, 
            annot=True,  # Show values in cells
            fmt=".1f",   # Format for annotations
            cmap="YlGnBu",
            linewidths=0.5,
            yticklabels=genes,
            xticklabels=data_cols,
            cbar_kws={"label": "Insertion Count"}
        )
        # Rotate column labels for better readability
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
    else:
        # Fallback to matplotlib
        plt.imshow(heatmap_data, aspect="auto", cmap="YlGnBu")
        plt.colorbar(label="Insertion Count")
        
        # Add labels
        plt.xticks(range(len(data_cols)), data_cols, rotation=45, ha="right")
        plt.yticks(range(len(genes)), genes)
        
        # Add values to cells
        for i in range(len(genes)):
            for j in range(len(data_cols)):
                plt.text(j, i, f"{heatmap_data[i, j]:.1f}", 
                         ha="center", va="center", 
                         color="black" if heatmap_data[i, j] < np.max(heatmap_data)/2 else "white")
    
    plt.title(title)
    plt.tight_layout()
    
    # Remove extension for base path
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created heatmap: {output_file} (PNG and PDF)")

def create_chromosome_heatmap(df, data_cols, output_file, title):
    """
    Create a chromosome heatmap visualization
    
    Args:
        df (pandas.DataFrame): Chromosome data frame
        data_cols (list): Data columns to include
        output_file (str): Output file path
        title (str): Chart title
    """
    # Extract data for heatmap
    heatmap_data = df[data_cols].values
    chromosomes = df["Chromosome"].values
    
    # Create figure with appropriate size
    plt.figure(figsize=(max(10, len(data_cols)), max(10, len(chromosomes)/2)))
    
    # Use seaborn if available for better visualization
    if HAS_SEABORN:
        ax = sns.heatmap(
            heatmap_data, 
            annot=True,  # Show values in cells
            fmt=".1f" if "normalized" in data_cols[0] else "g",  # Format for annotations
            cmap="YlGnBu",
            linewidths=0.5,
            yticklabels=chromosomes,
            xticklabels=data_cols,
            cbar_kws={"label": "Insertion Count" if "normalized" not in data_cols[0] else "% of Condition Total"}
        )
        # Rotate column labels for better readability
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
    else:
        # Fallback to matplotlib
        plt.imshow(heatmap_data, aspect="auto", cmap="YlGnBu")
        plt.colorbar(label="Insertion Count" if "normalized" not in data_cols[0] else "% of Condition Total")
        
        # Add labels
        plt.xticks(range(len(data_cols)), data_cols, rotation=45, ha="right")
        plt.yticks(range(len(chromosomes)), chromosomes)
        
        # Add values to cells
        for i in range(len(chromosomes)):
            for j in range(len(data_cols)):
                if "normalized" in data_cols[0]:
                    val_text = f"{heatmap_data[i, j]:.1f}"
                else:
                    val_text = f"{int(heatmap_data[i, j])}"
                
                plt.text(j, i, val_text, 
                         ha="center", va="center", 
                         color="black" if heatmap_data[i, j] < np.max(heatmap_data)/2 else "white")
    
    plt.title(title)
    plt.tight_layout()
    
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created chromosome heatmap: {output_file} (PNG and PDF)")

def create_condition_totals_chart(df, condition_cols, output_file):
    """
    Create a bar chart of total insertions per condition
    
    Args:
        df (pandas.DataFrame): Data frame
        condition_cols (list): Condition columns
        output_file (str): Output file path
    """
    # Calculate totals
    totals = [df[col].sum() for col in condition_cols]
    
    plt.figure(figsize=(max(8, len(condition_cols)/2), 6))
    
    # Create bar chart
    bars = plt.bar(condition_cols, totals, color="skyblue")
    
    # Add labels on top of bars
    for bar, total in zip(bars, totals):
        plt.text(
            bar.get_x() + bar.get_width()/2,
            bar.get_height() + (max(totals) * 0.01),
            f'{int(total):,}',
            ha='center', va='bottom',
            rotation=0
        )
    
    plt.title("Total Insertions by Condition")
    plt.ylabel("Total Insertion Count")
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created condition totals chart: {output_file} (PNG and PDF)")

def create_insertion_distribution(df, output_file):
    """
    Create a histogram of insertion counts per gene
    
    Args:
        df (pandas.DataFrame): Data frame
        output_file (str): Output file path
    """
    # Use Total column if available, otherwise sum the conditions
    if "Total" in df.columns:
        counts = df["Total"]
    else:
        metadata_cols = ["Gene", "GeneID"]
        data_cols = [col for col in df.columns if col not in metadata_cols]
        counts = df[data_cols].sum(axis=1)
    
    plt.figure(figsize=(10, 6))
    
    # Create histogram with log scale for better visualization
    if HAS_SEABORN:
        sns.histplot(counts, bins=30, kde=True)
    else:
        plt.hist(counts, bins=30, alpha=0.7, color="skyblue")
        
    plt.title("Distribution of Insertion Counts per Gene")
    plt.xlabel("Insertions per Gene")
    plt.ylabel("Number of Genes")
    plt.grid(linestyle="--", alpha=0.7)
    
    # Add log scale for better visualization if counts have a wide range
    if max(counts) / (min(counts) + 1) > 100:
        plt.xscale("log")
        plt.xlabel("Insertions per Gene (log scale)")
    
    plt.tight_layout()
    
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created insertion distribution chart: {output_file} (PNG and PDF)")

def create_gene_comparison_chart(df, data_cols, output_file, top_n=10, normalized=False):
    """
    Create a grouped bar chart comparing top genes across conditions
    
    Args:
        df (pandas.DataFrame): Data frame
        data_cols (list): Data columns to include
        output_file (str): Output file path
        top_n (int): Number of top genes to include
        normalized (bool): Whether the data is normalized
    """
    # Sort and select top genes
    if "Total" in df.columns:
        top_df = df.sort_values("Total", ascending=False).head(top_n)
    else:
        # If no Total column, sum the data columns
        sum_col = [col for col in data_cols if "normalized" not in col] if normalized else data_cols
        df["sum"] = df[sum_col].sum(axis=1)
        top_df = df.sort_values("sum", ascending=False).head(top_n)
    
    # Get gene names for legend
    genes = top_df["Gene"].values
    
    # Set up the figure with appropriate size
    plt.figure(figsize=(max(10, len(data_cols)*1.5), max(6, top_n*0.5)))
    
    # Calculate bar positions
    bar_width = 0.8 / top_n
    condition_positions = np.arange(len(data_cols))
    
    # Plot each gene as a grouped bar
    for i, gene in enumerate(genes):
        gene_data = top_df.loc[top_df["Gene"] == gene, data_cols].values[0]
        positions = condition_positions + (i - top_n/2 + 0.5) * bar_width
        
        plt.bar(
            positions, 
            gene_data, 
            width=bar_width, 
            label=gene
        )
    
    # Add labels and styling
    chart_title = "Top Genes Comparison (Normalized)" if normalized else "Top Genes Comparison"
    y_label = "% of Condition Total" if normalized else "Insertion Count"
    
    plt.title(chart_title)
    plt.ylabel(y_label)
    plt.xticks(condition_positions, data_cols, rotation=45, ha="right")
    plt.legend(title="Genes", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created gene comparison chart: {output_file} (PNG and PDF)")

def create_chromosome_profile_chart(df, data_cols, output_file, normalized=False):
    """
    Create line charts showing the insertion profile across chromosomes for each condition
    
    Args:
        df (pandas.DataFrame): Chromosome data frame
        data_cols (list): Data columns to include
        output_file (str): Output file path
        normalized (bool): Whether the data is normalized
    """
    # Set up the figure with appropriate size
    plt.figure(figsize=(12, 8))
    
    # Get chromosomes as x-axis values
    chromosomes = df["Chromosome"].values
    x = np.arange(len(chromosomes))
    
    # Plot each condition as a separate line
    for col in data_cols:
        plt.plot(x, df[col].values, marker='o', label=col.replace('_normalized', ''))
    
    # Add labels and styling
    chart_title = "Chromosome Insertion Profiles (Normalized)" if normalized else "Chromosome Insertion Profiles"
    y_label = "% of Condition Total" if normalized else "Insertion Count"
    
    plt.title(chart_title)
    plt.xlabel("Chromosome")
    plt.ylabel(y_label)
    plt.xticks(x, chromosomes, rotation=90)
    plt.legend(title="Conditions")
    plt.grid(linestyle="--", alpha=0.7)
    plt.tight_layout()
    
    # If normalized, set y-axis limit to better display variations
    if normalized:
        plt.ylim(0, min(25, df[data_cols].max().max() * 1.2))
    
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created chromosome profile chart: {output_file} (PNG and PDF)")

def create_chromosome_group_comparison(raw_df, norm_df, condition_cols, output_file):
    """
    Create a comparison of main chromosome groups across conditions
    
    Args:
        raw_df (pandas.DataFrame): Raw chromosome data
        norm_df (pandas.DataFrame): Normalized chromosome data
        condition_cols (list): Condition columns
        output_file (str): Output file path
    """
    # Define chromosome groups
    chrom_groups = {
        "Large (1-5)": ['1', '2', '3', '4', '5'],
        "Medium (6-12)": ['6', '7', '8', '9', '10', '11', '12'],
        "Small (13-22)": ['13', '14', '15', '16', '17', '18', '19', '20', '21', '22'],
        "Sex (X,Y)": ['X', 'Y']
    }
    
    # Create a dataframe for group statistics
    group_data = []
    
    for group_name, chromosomes in chrom_groups.items():
        group_row = {'Group': group_name}
        
        # Calculate raw counts for each condition
        for col in condition_cols:
            filtered_df = raw_df[raw_df['Chromosome'].isin(chromosomes)]
            group_row[col] = filtered_df[col].sum()
            
            # Add normalized percentage
            norm_col = f"{col}_normalized"
            filtered_norm_df = norm_df[norm_df['Chromosome'].isin(chromosomes)]
            group_row[norm_col] = filtered_norm_df[norm_col].sum()
        
        group_data.append(group_row)
    
    group_df = pd.DataFrame(group_data)
    
    # Create a figure with 2 subplots (raw and normalized)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
    
    # Plot raw counts
    group_df.plot(x='Group', y=condition_cols, kind='bar', ax=ax1)
    ax1.set_title("Raw Insertion Counts by Chromosome Group")
    ax1.set_ylabel("Insertion Count")
    ax1.legend(title="Conditions")
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Plot normalized counts
    normalized_cols = [f"{col}_normalized" for col in condition_cols]
    group_df.plot(x='Group', y=normalized_cols, kind='bar', ax=ax2)
    ax2.set_title("Normalized Insertion Distribution by Chromosome Group")
    ax2.set_ylabel("% of Condition Total")
    # Clean up legend labels
    ax2.legend(title="Conditions", labels=[col.replace('_normalized', '') for col in normalized_cols])
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    output_base = os.path.splitext(output_file)[0]
    save_figure_formats(plt, output_base)
    plt.close()
    
    log(f"Created chromosome group comparison: {output_file} (PNG and PDF)")

def save_normalized_matrix(norm_df, output_dir):
    """
    Save the normalized matrix to files
    
    Args:
        norm_df (pandas.DataFrame): Normalized matrix
        output_dir (str): Output directory
    """
    # Save full normalized matrix
    full_output = os.path.join(output_dir, "normalized_matrix_full.txt")
    norm_df.to_csv(full_output, sep='\t', index=False)
    log(f"Saved full normalized matrix to {full_output}")
    
    # Save a clean version with only original columns and normalized columns
    metadata_cols = ["Gene", "GeneID"]
    original_cols = [col for col in norm_df.columns if not col.endswith("_normalized") and col not in ["Total", "sum"]]
    normalized_cols = [col for col in norm_df.columns if col.endswith("_normalized")]
    
    # Combine metadata, original values, and normalized values
    clean_df = norm_df[metadata_cols + original_cols + normalized_cols]
    
    clean_output = os.path.join(output_dir, "normalized_matrix_clean.txt")
    clean_df.to_csv(clean_output, sep='\t', index=False)
    log(f"Saved clean normalized matrix to {clean_output}")
    
    # Create a split version for easier comparison
    condition_cols = [col for col in original_cols if col not in metadata_cols]
    
    # Write split matrix with side-by-side raw and normalized values
    split_output = os.path.join(output_dir, "normalized_matrix_split.txt")
    
    with open(split_output, 'w') as f:
        # Write header
        f.write("Gene\tGeneID")
        for col in condition_cols:
            f.write(f"\t{col}\t{col}_normalized")
        f.write("\n")
        
        # Write data
        for _, row in norm_df.iterrows():
            f.write(f"{row['Gene']}\t{row['GeneID']}")
            for col in condition_cols:
                norm_col = f"{col}_normalized"
                f.write(f"\t{row[col]}\t{row[norm_col]:.2f}")
            f.write("\n")
    
    log(f"Saved split normalized matrix to {split_output}")

def save_normalized_chromosome_matrix(norm_df, output_dir):
    """
    Save the normalized chromosome matrix to files
    
    Args:
        norm_df (pandas.DataFrame): Normalized chromosome matrix
        output_dir (str): Output directory
    """
    # Save full normalized matrix
    full_output = os.path.join(output_dir, "normalized_chromosome_matrix_full.txt")
    norm_df.to_csv(full_output, sep='\t', index=False)
    log(f"Saved full normalized chromosome matrix to {full_output}")
    
    # Save a clean version with only original columns and normalized columns
    original_cols = [col for col in norm_df.columns if not col.endswith("_normalized") and col != "sum"]
    normalized_cols = [col for col in norm_df.columns if col.endswith("_normalized")]
    
    # Combine original values and normalized values
    clean_df = norm_df[original_cols + normalized_cols]
    
    clean_output = os.path.join(output_dir, "normalized_chromosome_matrix_clean.txt")
    clean_df.to_csv(clean_output, sep='\t', index=False)
    log(f"Saved clean normalized chromosome matrix to {clean_output}")
    
    # Create a split version for easier comparison
    condition_cols = [col for col in original_cols if col != 'Chromosome']
    
    # Write split matrix with side-by-side raw and normalized values
    split_output = os.path.join(output_dir, "normalized_chromosome_matrix_split.txt")
    
    with open(split_output, 'w') as f:
        # Write header
        f.write("Chromosome")
        for col in condition_cols:
            f.write(f"\t{col}\t{col}_normalized")
        f.write("\n")
        
        # Write data
        for _, row in norm_df.iterrows():
            f.write(f"{row['Chromosome']}")
            for col in condition_cols:
                norm_col = f"{col}_normalized"
                f.write(f"\t{row[col]}\t{row[norm_col]:.2f}")
            f.write("\n")
    
    log(f"Saved split normalized chromosome matrix to {split_output}")

def create_significant_genes_excel(df, output_dir, top_n=50):
    """
    Create Excel file with significant gene information
    
    Args:
        df (pandas.DataFrame): Gene matrix data
        output_dir (str): Output directory for the Excel file
        top_n (int): Number of top genes to include (by Total insertions)
        
    Returns:
        str: Path to the created Excel file
    """
    if not HAS_XLSXWRITER:
        log("WARNING: xlsxwriter not available, cannot create Excel file")
        return None
    
    # Excel file path
    excel_file = os.path.join(output_dir, "significant_genes.xlsx")
    
    # Sort genes by Total insertions and take top N
    if "Total" in df.columns:
        top_genes = df.sort_values("Total", ascending=False).head(top_n)
    else:
        # If no Total column, sum the condition columns
        metadata_cols = ["Gene", "GeneID"]
        data_cols = [col for col in df.columns if col not in metadata_cols]
        df["Total"] = df[data_cols].sum(axis=1)
        top_genes = df.sort_values("Total", ascending=False).head(top_n)
    
    # Calculate normalized score (percentile within top genes)
    # This serves as a substitute for corrected p-value
    total_insertions = top_genes["Total"].sum()
    top_genes["Normalized_Score"] = (top_genes["Total"] / total_insertions) * 100
    
    # Create Excel workbook
    try:
        workbook = xlsxwriter.Workbook(excel_file)
        worksheet = workbook.add_worksheet("Significant Genes")
        
        # Add header formatting
        header_format = workbook.add_format({
            'bold': True,
            'fg_color': '#D7E4BC',
            'border': 1,
            'align': 'center',
        })
        
        # Add number formatting
        number_format = workbook.add_format({'num_format': '#,##0'})
        percent_format = workbook.add_format({'num_format': '0.00%'})
        
        # Add header row
        headers = ["Gene Name", "% of Total Insertions", "Total Insertions"]
        for col, header in enumerate(headers):
            worksheet.write(0, col, header, header_format)
        
        # Add gene data
        for row, (_, gene_data) in enumerate(top_genes.iterrows(), 1):
            worksheet.write(row, 0, gene_data["Gene"])
            # Using normalized score as a substitute for p-value
            worksheet.write(row, 1, gene_data["Normalized_Score"] / 100, percent_format)
            worksheet.write(row, 2, gene_data["Total"], number_format)
        
        # Adjust column widths
        worksheet.set_column(0, 0, 25)  # Gene Name
        worksheet.set_column(1, 1, 18)  # Corrected P-value
        worksheet.set_column(2, 2, 15)  # Total Insertions
        
        # Close workbook
        workbook.close()
        
        log(f"Created Excel file with {top_n} significant genes: {excel_file}")
        return excel_file
    
    except Exception as e:
        log(f"ERROR creating Excel file: {str(e)}")
        return None

def create_treatment_specific_excel(df, output_dir, top_n=50, heatmap_genes=20):
    """
    Create separate Excel files for each cordycepin condition
    
    Args:
        df (pandas.DataFrame): Gene matrix data
        output_dir (str): Output directory for the Excel files
        top_n (int): Number of top genes to include in Excel files
        heatmap_genes (int): Number of top genes to include in heatmaps
        
    Returns:
        list: Paths to the created Excel files
    """
    if not HAS_XLSXWRITER:
        log("WARNING: xlsxwriter not available, cannot create Excel files")
        return []
    
    # Define individual condition columns
    conditions = {
        "hyPB_cordycepin_150": "hyPB_Cordycepin_150uM",
        "hyPB_cordycepin_200": "hyPB_Cordycepin_200uM",
        "suPB_cordycepin_150": "suPB_Cordycepin_150uM",
        "suPB_cordycepin_200": "suPB_Cordycepin_200uM"
    }
    
    # Define control columns for each transposase type
    control_columns = {
        "hyPB": "hyPB_Control_DMSO",
        "suPB": "suPB_Control_DMSO"
    }
    
    excel_files = []
    
    # Process each condition separately
    for condition_name, column_name in conditions.items():
        # Check if this column exists
        if column_name not in df.columns:
            log(f"WARNING: Column '{column_name}' not found in data, skipping")
            continue
            
        log(f"Processing {condition_name} condition...")
        
        # Determine which control to use based on condition name
        transposase_type = "hyPB" if condition_name.startswith("hyPB") else "suPB"
        control_column = control_columns[transposase_type]
        
        # Check if control column exists
        if control_column not in df.columns:
            log(f"WARNING: Control column '{control_column}' not found in data, will exclude control data")
            include_control = False
        else:
            include_control = True
        
        # Create a copy to avoid modifying the original
        condition_df = df.copy()
        
        # Set the condition column as the total
        condition_df["Condition_Total"] = condition_df[column_name]
        
        # Calculate total insertions for the ENTIRE condition, not just top genes
        # This ensures consistency with the heatmap normalization
        total_condition_insertions = condition_df[column_name].sum()
        
        # Calculate normalized score for all genes first
        if total_condition_insertions > 0:
            condition_df["Normalized_Score"] = (condition_df[column_name] / total_condition_insertions) * 100
        else:
            condition_df["Normalized_Score"] = 0
            
        # Add control data normalization if control exists
        if include_control:
            # Store raw control values
            condition_df["Control_Total"] = condition_df[control_column]
            
            # Calculate normalized control data
            total_control_insertions = condition_df[control_column].sum()
            if total_control_insertions > 0:
                condition_df["Control_Normalized"] = (condition_df[control_column] / total_control_insertions) * 100
            else:
                condition_df["Control_Normalized"] = 0
        
        # Sort and filter AFTER normalization
        condition_genes = condition_df.sort_values("Condition_Total", ascending=False).head(top_n)
            
        # Create Excel file
        excel_file = os.path.join(output_dir, f"{condition_name}_significant_genes.xlsx")
        try:
            workbook = xlsxwriter.Workbook(excel_file)
            worksheet = workbook.add_worksheet(f"{condition_name}")
            
            # Add header formatting
            header_format = workbook.add_format({
                'bold': True,
                'fg_color': '#D7E4BC',
                'border': 1,
                'align': 'center',
            })
            
            # Add number formatting
            number_format = workbook.add_format({'num_format': '#,##0'})
            percent_format = workbook.add_format({'num_format': '0.00%'})
            
            # Prepare headers based on whether control data is included
            if include_control:
                headers = [
                    "Gene Name", 
                    "% of Total Insertions", 
                    "Insertions",
                    f"{control_column} Insertions",
                    f"{control_column} %"
                ]
            else:
                headers = ["Gene Name", "% of Total Insertions", "Insertions"]
                
            # Add header row
            for col, header in enumerate(headers):
                worksheet.write(0, col, header, header_format)
            
            # Add gene data
            for row, (_, gene_data) in enumerate(condition_genes.iterrows(), 1):
                worksheet.write(row, 0, gene_data["Gene"])
                # Using normalized score as a percentage
                worksheet.write(row, 1, gene_data["Normalized_Score"] / 100, percent_format)
                worksheet.write(row, 2, gene_data["Condition_Total"], number_format)
                
                # Add control data if available
                if include_control:
                    worksheet.write(row, 3, gene_data["Control_Total"], number_format)
                    worksheet.write(row, 4, gene_data["Control_Normalized"] / 100, percent_format)
            
            # Adjust column widths
            worksheet.set_column(0, 0, 25)  # Gene Name
            worksheet.set_column(1, 1, 18)  # % of Total Insertions
            worksheet.set_column(2, 2, 15)  # Insertions
            if include_control:
                worksheet.set_column(3, 3, 25)  # Control Insertions
                worksheet.set_column(4, 4, 18)  # Control %
            
            # Close workbook
            workbook.close()
            
            log(f"Created {condition_name} Excel file with {'control data' if include_control else 'no control data'}: {excel_file}")
            excel_files.append(excel_file)
            
            # Create heatmap if control data is available
            if include_control:
                create_excel_heatmap(condition_genes, condition_name, control_column, output_dir, heatmap_genes)
            
        except Exception as e:
            log(f"ERROR creating {condition_name} Excel file: {str(e)}")
    
    if excel_files:
        log(f"Created {len(excel_files)} condition-specific Excel files")
    else:
        log("WARNING: No condition-specific Excel files were created")
        
    return excel_files

def create_excel_heatmap(gene_df, condition_name, control_column, output_dir, top_n=20):
    """
    Create a heatmap for treatment vs. control data from Excel files
    
    Args:
        gene_df (pandas.DataFrame): Gene data with treatment and control values
        condition_name (str): Name of the treatment condition
        control_column (str): Name of the control column
        output_dir (str): Directory to save the heatmap
        top_n (int): Number of top genes to include in the heatmap
    """
    # Check if we have both treatment and control data
    if "Normalized_Score" not in gene_df.columns or "Control_Normalized" not in gene_df.columns:
        log(f"WARNING: Missing normalized data for heatmap creation, skipping {condition_name} heatmap")
        return None
    
    try:
        # Create visualizations directory if it doesn't exist
        viz_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)
        
        # File path for the heatmap
        heatmap_file = os.path.join(viz_dir, f"{condition_name}_vs_control_heatmap.png")
        
        # Take top N genes for better visualization (or fewer if we have less)
        genes_to_show = min(top_n, len(gene_df))
        top_genes = gene_df.head(genes_to_show)
        
        # Extract data for heatmap
        gene_names = top_genes["Gene"].values
        treatment_values = top_genes["Normalized_Score"].values
        control_values = top_genes["Control_Normalized"].values
        
        # Create a 2D array for the heatmap
        heatmap_data = np.column_stack((treatment_values, control_values))
        
        # Create figure with appropriate size
        plt.figure(figsize=(8, max(8, len(gene_names)/1.5)))
        
        # Column labels for the heatmap
        column_labels = [condition_name, control_column]
        
        # Use seaborn if available for better visualization
        if HAS_SEABORN:
            ax = sns.heatmap(
                heatmap_data, 
                annot=True,  # Show values in cells
                fmt=".2f",   # Format for annotations
                cmap="YlGnBu",
                linewidths=0.5,
                yticklabels=gene_names,
                xticklabels=column_labels,
                cbar_kws={"label": "% of Total Insertions"}
            )
            # Rotate column labels for better readability
            plt.xticks(rotation=45, ha="right")
            plt.yticks(rotation=0)
        else:
            # Fallback to matplotlib
            plt.imshow(heatmap_data, aspect="auto", cmap="YlGnBu")
            plt.colorbar(label="% of Total Insertions")
            
            # Add labels
            plt.xticks(range(len(column_labels)), column_labels, rotation=45, ha="right")
            plt.yticks(range(len(gene_names)), gene_names)
            
            # Add values to cells
            for i in range(len(gene_names)):
                for j in range(len(column_labels)):
                    plt.text(j, i, f"{heatmap_data[i, j]:.2f}", 
                             ha="center", va="center", 
                             color="black" if heatmap_data[i, j] < np.max(heatmap_data)/2 else "white")
        
        plt.title(f"{condition_name} vs {control_column} Comparison")
        plt.tight_layout()
        
        heatmap_file = os.path.join(viz_dir, f"{condition_name}_vs_control_heatmap.png")
        output_base = os.path.splitext(heatmap_file)[0]
        save_figure_formats(plt, output_base)
        plt.close()
        
        log(f"Created heatmap for {condition_name} vs control: {heatmap_file} (PNG and PDF)")
        return heatmap_file
    
    except Exception as e:
        log(f"ERROR creating heatmap for {condition_name}: {str(e)}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Creates normalized and improved visualizations for condition comparison matrices."
    )
    parser.add_argument("--matrix", required=True,
        help="Input gene-condition matrix file (TSV format).")
    parser.add_argument("--output_dir", required=True,
        help="Output directory for normalized matrix and visualizations.")
    parser.add_argument("--top_genes", type=int, default=50,
        help="Number of top genes to include in Excel files (default: 50).")
    parser.add_argument("--heatmap_genes", type=int, default=20,
        help="Number of top genes to include in heatmaps (default: 20).")
    parser.add_argument("--chrom_matrix", default=None,
        help="Optional chromosome-condition matrix file (TSV format).")
    parser.add_argument("--excel", action="store_true",
        help="Generate Excel file with significant genes.")

    args = parser.parse_args()
    start_time = time.time()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process gene matrix
    log("Processing gene matrix...")
    df = load_matrix(args.matrix)
    if df is not None:
        # Create normalized matrix
        norm_df = create_normalized_matrix(df)
        if norm_df is not None:
            # Save normalized matrices
            save_normalized_matrix(norm_df, args.output_dir)
            
            # Create visualizations
            viz_dir = os.path.join(args.output_dir, "visualizations")
            create_visualizations(df, norm_df, viz_dir, top_n=args.top_genes, heatmap_genes=args.heatmap_genes)
            
            # Create Excel file with significant genes
            if args.excel and HAS_XLSXWRITER:
                # Original Excel with all conditions
                create_significant_genes_excel(df, args.output_dir, args.top_genes)
                # Treatment-specific Excel files
                create_treatment_specific_excel(df, args.output_dir, args.top_genes, args.heatmap_genes)
            elif args.excel and not HAS_XLSXWRITER:
                log("WARNING: Excel generation requested but xlsxwriter library not found. Please install with 'pip install xlsxwriter'")
        else:
            log("ERROR: Failed to create normalized gene matrix. Skipping visualizations.")
    else:
        log("ERROR: Failed to load gene matrix. Skipping gene matrix analysis.")
    
    # Process chromosome matrix if provided
    if args.chrom_matrix:
        log("Processing chromosome matrix...")
        chrom_df = load_chromosome_matrix(args.chrom_matrix)
        if chrom_df is not None:
            # Create normalized chromosome matrix
            norm_chrom_df = create_normalized_chromosome_matrix(chrom_df)
            if norm_chrom_df is not None:
                # Save normalized chromosome matrices
                save_normalized_chromosome_matrix(norm_chrom_df, args.output_dir)
                
                # Create chromosome visualizations
                chrom_viz_dir = os.path.join(args.output_dir, "chromosome_visualizations")
                create_chromosome_visualizations(chrom_df, norm_chrom_df, chrom_viz_dir)
            else:
                log("ERROR: Failed to create normalized chromosome matrix. Skipping visualizations.")
        else:
            log("ERROR: Failed to load chromosome matrix. Skipping chromosome analysis.")
    
    # Report completion
    elapsed = time.time() - start_time
    log(f"Matrix normalization and visualization complete in {elapsed:.1f} seconds")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"ERROR: {str(e)}")
        sys.exit(1)
