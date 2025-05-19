#!/usr/bin/env python3
"""
Normalize_Matrix_Visualization.py
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

Usage Example:
  python Normalize_Matrix_Visualization.py \
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

def create_visualizations(df, norm_df, output_dir, top_n=50):
    """
    Create matrix visualizations
    
    Args:
        df (pandas.DataFrame): Original matrix
        norm_df (pandas.DataFrame): Normalized matrix
        output_dir (str): Output directory for visualizations
        top_n (int): Number of top genes to include in heatmaps
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
            top_n=top_n
        )
    except Exception as e:
        log(f"ERROR creating raw counts heatmap: {str(e)}")
    
    # 2. Normalized counts heatmap for top genes
    try:
        create_heatmap(
            norm_df, normalized_cols, 
            os.path.join(output_dir, "top_genes_heatmap_normalized.png"),
            "Normalized Insertion Counts (% of condition total)",
            top_n=top_n
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
    
    # 5. Gene counts across conditions (for top 10 genes)
    try:
        create_gene_comparison_chart(
            df, condition_cols,
            os.path.join(output_dir, "top_genes_comparison.png"),
            top_n=10
        )
    except Exception as e:
        log(f"ERROR creating gene comparison chart: {str(e)}")
    
    # 6. Normalized gene counts (for same top 10 genes)
    try:
        create_gene_comparison_chart(
            norm_df, normalized_cols,
            os.path.join(output_dir, "top_genes_comparison_normalized.png"),
            top_n=10,
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
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created heatmap: {output_file}")

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
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created chromosome heatmap: {output_file}")

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
    
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created condition totals chart: {output_file}")

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
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created insertion distribution chart: {output_file}")

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
    
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created gene comparison chart: {output_file}")

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
    
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created chromosome profile chart: {output_file}")

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
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    log(f"Created chromosome group comparison: {output_file}")

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

def main():
    parser = argparse.ArgumentParser(
        description="Creates normalized and improved visualizations for condition comparison matrices."
    )
    parser.add_argument("--matrix", required=True,
        help="Input gene-condition matrix file (TSV format).")
    parser.add_argument("--output_dir", required=True,
        help="Output directory for normalized matrix and visualizations.")
    parser.add_argument("--top_genes", type=int, default=50,
        help="Number of top genes to include in visualizations (default: 50).")
    parser.add_argument("--chrom_matrix", default=None,
        help="Optional chromosome-condition matrix file (TSV format).")

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
            create_visualizations(df, norm_df, viz_dir, top_n=args.top_genes)
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
