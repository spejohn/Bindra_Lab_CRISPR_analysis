"""
Quality control plotting functions for CRISPR screening analysis.
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any, Set
from matplotlib.figure import Figure

# Import configurations
from analysis_pipeline.core.config import DEFAULT_ESSENTIAL_GENES, DEFAULT_FDR_THRESHOLD
from analysis_pipeline.core.file_handling import ensure_output_dir


def plot_count_distribution(
    count_table: str,
    output_dir: str,
    log_scale: bool = True,
    min_cutoff: int = 0,
    figsize: Tuple[int, int] = (12, 8),
    output_prefix: str = "count_distribution"
) -> Dict[str, str]:
    """
    Plot distribution of guide counts across samples.
    
    Args:
        count_table: Path to count table file
        output_dir: Directory to save plot files
        log_scale: Whether to use log scale for counts
        min_cutoff: Minimum count value to include
        figsize: Figure size as (width, height)
        output_prefix: Prefix for output files
        
    Returns:
        Dictionary of output file paths
    """
    logging.info(f"Plotting count distribution from {count_table}")
    
    # Ensure output directory exists
    ensure_output_dir(output_dir)
    
    # Read count table
    try:
        df = pd.read_csv(count_table, sep='\t')
        logging.info(f"Read count table with {len(df)} rows")
    except Exception as e:
        error_msg = f"Error reading count table {count_table}: {e}"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Identify guide and count columns
    guide_column = None
    for col in ["sgRNA", "sgrna", "guide"]:
        if col in df.columns:
            guide_column = col
            break
            
    if guide_column is None:
        error_msg = "Could not identify guide column in count table"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Identify count columns (exclude guide, gene columns)
    exclude_cols = [guide_column]
    for col in ["Gene", "gene"]:
        if col in df.columns:
            exclude_cols.append(col)
    
    count_columns = [col for col in df.columns if col not in exclude_cols]
    if not count_columns:
        error_msg = "No count columns found in count table"
        logging.error(error_msg)
        return {"error": error_msg}
    
    logging.info(f"Identified {len(count_columns)} count columns: {', '.join(count_columns)}")
    
    # Filter counts by minimum cutoff
    if min_cutoff > 0:
        logging.info(f"Filtering counts with minimum value of {min_cutoff}")
        
    # Create output files dictionary
    output_files = {}
    
    # Plot 1: Count distributions as boxplot
    fig, ax = plt.subplots(figsize=figsize)
    plot_data = []
    
    for col in count_columns:
        counts = df[col].values
        if min_cutoff > 0:
            counts = counts[counts >= min_cutoff]
        
        if log_scale:
            # Handle zeros when using log scale
            counts = np.log10(counts + 1)
        
        plot_data.append(counts)
    
    ax.boxplot(plot_data, labels=count_columns)
    ax.set_title("Guide Count Distribution by Sample")
    ax.set_ylabel("Count" if not log_scale else "Log10(Count + 1)")
    ax.set_xticklabels(count_columns, rotation=90)
    
    plt.tight_layout()
    boxplot_file = os.path.join(output_dir, f"{output_prefix}_boxplot.png")
    fig.savefig(boxplot_file, dpi=300)
    plt.close(fig)
    
    output_files["boxplot"] = boxplot_file
    logging.info(f"Saved count distribution boxplot to {boxplot_file}")
    
    # Plot 2: Count distributions as violin plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create a DataFrame in long format for seaborn
    long_df = pd.DataFrame(columns=["Sample", "Count"])
    
    for col in count_columns:
        counts = df[col].values
        if min_cutoff > 0:
            counts = counts[counts >= min_cutoff]
        
        if log_scale:
            # Handle zeros when using log scale
            counts = np.log10(counts + 1)
        
        sample_df = pd.DataFrame({
            "Sample": [col] * len(counts),
            "Count": counts
        })
        
        long_df = pd.concat([long_df, sample_df], ignore_index=True)
    
    sns.violinplot(x="Sample", y="Count", data=long_df, ax=ax)
    ax.set_title("Guide Count Distribution by Sample")
    ax.set_ylabel("Count" if not log_scale else "Log10(Count + 1)")
    ax.set_xticklabels(count_columns, rotation=90)
    
    plt.tight_layout()
    violin_file = os.path.join(output_dir, f"{output_prefix}_violin.png")
    fig.savefig(violin_file, dpi=300)
    plt.close(fig)
    
    output_files["violin"] = violin_file
    logging.info(f"Saved count distribution violin plot to {violin_file}")
    
    # Plot 3: Histogram of counts for each sample
    fig, axes = plt.subplots(
        nrows=(len(count_columns) + 2) // 3, 
        ncols=3, 
        figsize=(15, 5 * ((len(count_columns) + 2) // 3)),
        sharex=True
    )
    axes = axes.flatten()
    
    for i, col in enumerate(count_columns):
        if i < len(axes):
            counts = df[col].values
            if min_cutoff > 0:
                counts = counts[counts >= min_cutoff]
            
            if log_scale:
                # Use log-scale bins
                bins = np.logspace(0, np.log10(max(counts) + 1), 50)
                counts_plot = counts
            else:
                bins = 50
                counts_plot = counts
            
            axes[i].hist(counts_plot, bins=bins)
            axes[i].set_title(col)
            
            if log_scale:
                axes[i].set_xscale('log')
            
            axes[i].set_xlabel("Count")
            axes[i].set_ylabel("Frequency")
    
    # Hide unused subplots
    for i in range(len(count_columns), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    hist_file = os.path.join(output_dir, f"{output_prefix}_histogram.png")
    fig.savefig(hist_file, dpi=300)
    plt.close(fig)
    
    output_files["histogram"] = hist_file
    logging.info(f"Saved count histograms to {hist_file}")
    
    return output_files


def plot_guide_correlations(
    count_table: str,
    output_dir: str,
    log_transform: bool = True,
    min_cutoff: int = 0,
    sample_subset: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (12, 10),
    output_prefix: str = "guide_correlation"
) -> Dict[str, str]:
    """
    Plot correlation of guide counts between samples.
    
    Args:
        count_table: Path to count table file
        output_dir: Directory to save plot files
        log_transform: Whether to log-transform counts
        min_cutoff: Minimum count to include
        sample_subset: Subset of samples to plot (if None, use all)
        figsize: Figure size as (width, height)
        output_prefix: Prefix for output files
        
    Returns:
        Dictionary of output file paths
    """
    logging.info(f"Plotting guide count correlations from {count_table}")
    
    # Ensure output directory exists
    ensure_output_dir(output_dir)
    
    # Read count table
    try:
        df = pd.read_csv(count_table, sep='\t')
        logging.info(f"Read count table with {len(df)} rows")
    except Exception as e:
        error_msg = f"Error reading count table {count_table}: {e}"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Identify guide and count columns
    guide_column = None
    for col in ["sgRNA", "sgrna", "guide"]:
        if col in df.columns:
            guide_column = col
            break
            
    if guide_column is None:
        error_msg = "Could not identify guide column in count table"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Identify count columns (exclude guide, gene columns)
    exclude_cols = [guide_column]
    for col in ["Gene", "gene"]:
        if col in df.columns:
            exclude_cols.append(col)
    
    count_columns = [col for col in df.columns if col not in exclude_cols]
    if not count_columns:
        error_msg = "No count columns found in count table"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Filter by sample subset if provided
    if sample_subset:
        count_columns = [col for col in count_columns if col in sample_subset]
        if not count_columns:
            error_msg = "No samples from subset found in count table"
            logging.error(error_msg)
            return {"error": error_msg}
    
    logging.info(f"Using {len(count_columns)} count columns for correlation analysis")
    
    # Create output files dictionary
    output_files = {}
    
    # Extract only count columns for correlation
    counts_df = df[count_columns].copy()
    
    # Apply log transformation if requested
    if log_transform:
        logging.info("Applying log10 transformation to counts")
        counts_df = np.log10(counts_df + 1)
    
    # Apply minimum cutoff
    if min_cutoff > 0:
        logging.info(f"Filtering counts with minimum value of {min_cutoff}")
        for col in count_columns:
            if log_transform:
                # Convert min_cutoff to log space
                log_cutoff = np.log10(min_cutoff + 1)
                counts_df = counts_df[counts_df[col] >= log_cutoff]
            else:
                counts_df = counts_df[counts_df[col] >= min_cutoff]
    
    # Calculate correlation matrix
    corr_matrix = counts_df.corr()
    
    # Plot 1: Correlation heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(
        corr_matrix, 
        annot=True, 
        cmap="coolwarm", 
        vmin=-1, 
        vmax=1, 
        ax=ax
    )
    
    ax.set_title("Guide Count Correlation Between Samples")
    
    plt.tight_layout()
    heatmap_file = os.path.join(output_dir, f"{output_prefix}_heatmap.png")
    fig.savefig(heatmap_file, dpi=300)
    plt.close(fig)
    
    output_files["heatmap"] = heatmap_file
    logging.info(f"Saved correlation heatmap to {heatmap_file}")
    
    # Plot 2: Scatter plots for pairs of samples
    # If there are many samples, limit to first few pairs
    if len(count_columns) > 5:
        logging.info(f"Too many samples ({len(count_columns)}) for pairwise plots, limiting to first 5")
        plot_columns = count_columns[:5]
    else:
        plot_columns = count_columns
    
    # Create pairwise scatter plots
    for i, col1 in enumerate(plot_columns):
        for j, col2 in enumerate(plot_columns):
            if i >= j:  # Skip duplicates and self-comparisons
                continue
                
            fig, ax = plt.subplots(figsize=(8, 8))
            
            # Get correlation coefficient
            corr = corr_matrix.loc[col1, col2]
            
            # Plot scatter with hexbin for dense data
            if len(df) > 5000:
                hb = ax.hexbin(
                    counts_df[col1], 
                    counts_df[col2], 
                    gridsize=50, 
                    cmap='viridis',
                    mincnt=1
                )
                plt.colorbar(hb, ax=ax, label='Count')
            else:
                ax.scatter(
                    counts_df[col1], 
                    counts_df[col2], 
                    alpha=0.5,
                    s=5
                )
            
            # Add correlation line
            m, b = np.polyfit(counts_df[col1], counts_df[col2], 1)
            ax.plot(counts_df[col1], m*counts_df[col1] + b, color='red', linestyle='--')
            
            # Add correlation coefficient to title
            ax.set_title(f"Correlation: {corr:.4f}")
            
            # Labels
            count_label = "Log10(Count + 1)" if log_transform else "Count"
            ax.set_xlabel(f"{col1} {count_label}")
            ax.set_ylabel(f"{col2} {count_label}")
            
            plt.tight_layout()
            scatter_file = os.path.join(output_dir, f"{output_prefix}_scatter_{col1}_vs_{col2}.png")
            fig.savefig(scatter_file, dpi=300)
            plt.close(fig)
            
            output_files[f"scatter_{col1}_vs_{col2}"] = scatter_file
            logging.info(f"Saved scatter plot for {col1} vs {col2} to {scatter_file}")
    
    # Save correlation matrix as CSV
    corr_file = os.path.join(output_dir, f"{output_prefix}_correlation_matrix.csv")
    corr_matrix.to_csv(corr_file)
    output_files["correlation_matrix"] = corr_file
    logging.info(f"Saved correlation matrix to {corr_file}")
    
    return output_files


def plot_mageck_results(
    gene_summary_file: str,
    output_dir: str,
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
    essential_genes: Optional[Dict[str, List[str]]] = None,
    output_prefix: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8)
) -> Dict[str, str]:
    """
    Plot MAGeCK analysis results.
    
    Args:
        gene_summary_file: Path to MAGeCK gene summary file
        output_dir: Directory to save plot files
        fdr_threshold: FDR threshold for significance
        essential_genes: Dictionary of essential gene lists
        output_prefix: Prefix for output files (default: derived from input)
        figsize: Figure size as (width, height)
        
    Returns:
        Dictionary of output file paths
    """
    logging.info(f"Plotting MAGeCK results from {gene_summary_file}")
    
    # Ensure output directory exists
    ensure_output_dir(output_dir)
    
    # Set default output prefix based on input file if not provided
    if output_prefix is None:
        output_prefix = Path(gene_summary_file).stem.replace(".gene_summary", "")
    
    # Read gene summary file
    try:
        df = pd.read_csv(gene_summary_file, sep='\t')
        logging.info(f"Read gene summary with {len(df)} rows")
    except Exception as e:
        error_msg = f"Error reading gene summary file {gene_summary_file}: {e}"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Check if required columns exist
    required_columns = ["id", "neg|lfc", "neg|fdr", "pos|lfc", "pos|fdr"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        error_msg = f"Gene summary file is missing required columns: {', '.join(missing_columns)}"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Use default essential genes if not provided
    if essential_genes is None:
        essential_genes = DEFAULT_ESSENTIAL_GENES
    
    # Create a set of all essential genes
    all_essential_genes = set()
    for gene_list in essential_genes.values():
        all_essential_genes.update(gene_list)
    
    logging.info(f"Using {len(all_essential_genes)} essential genes for QC")
    
    # Create output files dictionary
    output_files = {}
    
    # Prepare data for plotting
    df["neg_significant"] = df["neg|fdr"] < fdr_threshold
    df["pos_significant"] = df["pos|fdr"] < fdr_threshold
    df["essential"] = df["id"].isin(all_essential_genes)
    
    # Count significant genes
    neg_sig_count = df["neg_significant"].sum()
    pos_sig_count = df["pos_significant"].sum()
    logging.info(f"Found {neg_sig_count} negative significant genes and {pos_sig_count} positive significant genes")
    
    # Count essential genes
    essential_count = df["essential"].sum()
    neg_essential_sig = ((df["essential"] == True) & (df["neg_significant"] == True)).sum()
    pos_essential_sig = ((df["essential"] == True) & (df["pos_significant"] == True)).sum()
    
    logging.info(f"Found {essential_count} essential genes in the dataset")
    logging.info(f"  {neg_essential_sig} essential genes are significant in negative selection")
    logging.info(f"  {pos_essential_sig} essential genes are significant in positive selection")
    
    # Plot 1: Volcano plot for negative selection
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot non-significant genes
    ax.scatter(
        df.loc[~df["neg_significant"] & ~df["essential"], "neg|lfc"],
        -np.log10(df.loc[~df["neg_significant"] & ~df["essential"], "neg|fdr"]),
        alpha=0.5,
        s=5,
        color="gray",
        label="Non-significant"
    )
    
    # Plot significant non-essential genes
    ax.scatter(
        df.loc[df["neg_significant"] & ~df["essential"], "neg|lfc"],
        -np.log10(df.loc[df["neg_significant"] & ~df["essential"], "neg|fdr"]),
        alpha=0.7,
        s=10,
        color="blue",
        label="Significant"
    )
    
    # Plot essential genes
    ax.scatter(
        df.loc[df["essential"], "neg|lfc"],
        -np.log10(df.loc[df["essential"], "neg|fdr"]),
        alpha=0.7,
        s=15,
        color="red",
        label="Essential"
    )
    
    # Add threshold line
    ax.axhline(-np.log10(fdr_threshold), linestyle="--", color="black", alpha=0.5)
    
    # Labels
    ax.set_title(f"Negative Selection (Dropouts)\n{neg_sig_count} significant genes, {neg_essential_sig} essential")
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel("-log10(FDR)")
    ax.legend()
    
    plt.tight_layout()
    neg_volcano_file = os.path.join(output_dir, f"{output_prefix}_neg_volcano.png")
    fig.savefig(neg_volcano_file, dpi=300)
    plt.close(fig)
    
    output_files["neg_volcano"] = neg_volcano_file
    logging.info(f"Saved negative selection volcano plot to {neg_volcano_file}")
    
    # Plot 2: Volcano plot for positive selection
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot non-significant genes
    ax.scatter(
        df.loc[~df["pos_significant"] & ~df["essential"], "pos|lfc"],
        -np.log10(df.loc[~df["pos_significant"] & ~df["essential"], "pos|fdr"]),
        alpha=0.5,
        s=5,
        color="gray",
        label="Non-significant"
    )
    
    # Plot significant non-essential genes
    ax.scatter(
        df.loc[df["pos_significant"] & ~df["essential"], "pos|lfc"],
        -np.log10(df.loc[df["pos_significant"] & ~df["essential"], "pos|fdr"]),
        alpha=0.7,
        s=10,
        color="green",
        label="Significant"
    )
    
    # Plot essential genes
    ax.scatter(
        df.loc[df["essential"], "pos|lfc"],
        -np.log10(df.loc[df["essential"], "pos|fdr"]),
        alpha=0.7,
        s=15,
        color="red",
        label="Essential"
    )
    
    # Add threshold line
    ax.axhline(-np.log10(fdr_threshold), linestyle="--", color="black", alpha=0.5)
    
    # Labels
    ax.set_title(f"Positive Selection (Enrichment)\n{pos_sig_count} significant genes, {pos_essential_sig} essential")
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel("-log10(FDR)")
    ax.legend()
    
    plt.tight_layout()
    pos_volcano_file = os.path.join(output_dir, f"{output_prefix}_pos_volcano.png")
    fig.savefig(pos_volcano_file, dpi=300)
    plt.close(fig)
    
    output_files["pos_volcano"] = pos_volcano_file
    logging.info(f"Saved positive selection volcano plot to {pos_volcano_file}")
    
    # Plot 3: Combined negative and positive LFC
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot non-significant genes
    ax.scatter(
        df.loc[~df["neg_significant"] & ~df["pos_significant"] & ~df["essential"], "neg|lfc"],
        df.loc[~df["neg_significant"] & ~df["pos_significant"] & ~df["essential"], "pos|lfc"],
        alpha=0.5,
        s=5,
        color="gray",
        label="Non-significant"
    )
    
    # Plot negative significant genes
    ax.scatter(
        df.loc[df["neg_significant"] & ~df["pos_significant"] & ~df["essential"], "neg|lfc"],
        df.loc[df["neg_significant"] & ~df["pos_significant"] & ~df["essential"], "pos|lfc"],
        alpha=0.7,
        s=10,
        color="blue",
        label="Negative significant"
    )
    
    # Plot positive significant genes
    ax.scatter(
        df.loc[~df["neg_significant"] & df["pos_significant"] & ~df["essential"], "neg|lfc"],
        df.loc[~df["neg_significant"] & df["pos_significant"] & ~df["essential"], "pos|lfc"],
        alpha=0.7,
        s=10,
        color="green",
        label="Positive significant"
    )
    
    # Plot both significant genes
    ax.scatter(
        df.loc[df["neg_significant"] & df["pos_significant"] & ~df["essential"], "neg|lfc"],
        df.loc[df["neg_significant"] & df["pos_significant"] & ~df["essential"], "pos|lfc"],
        alpha=0.7,
        s=10,
        color="purple",
        label="Both significant"
    )
    
    # Plot essential genes
    ax.scatter(
        df.loc[df["essential"], "neg|lfc"],
        df.loc[df["essential"], "pos|lfc"],
        alpha=0.7,
        s=15,
        color="red",
        label="Essential"
    )
    
    # Add reference lines
    ax.axhline(0, linestyle="--", color="black", alpha=0.3)
    ax.axvline(0, linestyle="--", color="black", alpha=0.3)
    
    # Labels
    ax.set_title(f"Combined LFC Comparison\n{essential_count} essential genes")
    ax.set_xlabel("Negative Selection LFC")
    ax.set_ylabel("Positive Selection LFC")
    ax.legend()
    
    plt.tight_layout()
    combined_file = os.path.join(output_dir, f"{output_prefix}_combined_lfc.png")
    fig.savefig(combined_file, dpi=300)
    plt.close(fig)
    
    output_files["combined_lfc"] = combined_file
    logging.info(f"Saved combined LFC plot to {combined_file}")
    
    return output_files


def plot_essential_gene_enrichment(
    gene_summary_files: Dict[str, str],
    output_dir: str,
    essential_genes: Optional[Dict[str, List[str]]] = None,
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
    figsize: Tuple[int, int] = (12, 10),
    output_prefix: str = "essential_enrichment"
) -> Dict[str, str]:
    """
    Plot essential gene enrichment across multiple contrasts.
    
    Args:
        gene_summary_files: Dictionary mapping contrast names to gene summary files
        output_dir: Directory to save plot files
        essential_genes: Dictionary of essential gene lists
        fdr_threshold: FDR threshold for significance
        figsize: Figure size as (width, height)
        output_prefix: Prefix for output files
        
    Returns:
        Dictionary of output file paths
    """
    logging.info(f"Plotting essential gene enrichment for {len(gene_summary_files)} contrasts")
    
    # Use default essential genes if not provided
    if essential_genes is None:
        essential_genes = DEFAULT_ESSENTIAL_GENES
    
    # Ensure output directory exists
    ensure_output_dir(output_dir)
    
    # Create a set of all essential genes
    all_essential_genes = set()
    for gene_list in essential_genes.values():
        all_essential_genes.update(gene_list)
    
    # Create output files dictionary
    output_files = {}
    
    # Process each gene summary file
    results = {}
    
    for contrast, file_path in gene_summary_files.items():
        try:
            # Read gene summary file
            df = pd.read_csv(file_path, sep='\t')
            logging.info(f"Read gene summary for {contrast} with {len(df)} rows")
            
            # Check if required columns exist
            required_columns = ["id", "neg|lfc", "neg|fdr"]
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                logging.error(f"Gene summary file for {contrast} is missing required columns: {', '.join(missing_columns)}")
                continue
            
            # Mark essential genes
            df["essential"] = df["id"].isin(all_essential_genes)
            
            # Count significant genes
            total_genes = len(df)
            total_essential = df["essential"].sum()
            
            neg_sig = (df["neg|fdr"] < fdr_threshold).sum()
            neg_essential_sig = ((df["essential"] == True) & (df["neg|fdr"] < fdr_threshold)).sum()
            
            # Store results
            results[contrast] = {
                "total_genes": total_genes,
                "total_essential": total_essential,
                "neg_sig": neg_sig,
                "neg_essential_sig": neg_essential_sig,
                "neg_sig_percent": 100 * neg_sig / total_genes,
                "neg_essential_percent": 100 * neg_essential_sig / total_essential if total_essential > 0 else 0
            }
            
            logging.info(f"Contrast {contrast}: {neg_sig} significant genes, {neg_essential_sig} essential")
            
        except Exception as e:
            logging.error(f"Error processing gene summary file for {contrast}: {e}")
            continue
    
    if not results:
        error_msg = "No valid gene summary files processed"
        logging.error(error_msg)
        return {"error": error_msg}
    
    # Create a DataFrame from results
    results_df = pd.DataFrame.from_dict(results, orient='index')
    
    # Sort by essential gene enrichment
    results_df = results_df.sort_values("neg_essential_percent", ascending=False)
    
    # Plot 1: Essential gene enrichment bar chart
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot bars for each contrast
    bar_width = 0.35
    x = np.arange(len(results_df))
    
    # Plot percentage of significant genes
    ax.bar(
        x - bar_width/2, 
        results_df["neg_sig_percent"], 
        width=bar_width, 
        label="All Genes",
        color="gray"
    )
    
    # Plot percentage of significant essential genes
    ax.bar(
        x + bar_width/2, 
        results_df["neg_essential_percent"], 
        width=bar_width, 
        label="Essential Genes",
        color="red"
    )
    
    # Add labels
    ax.set_title("Essential Gene Enrichment in Negative Selection")
    ax.set_xlabel("Contrast")
    ax.set_ylabel("Percent Significant (FDR < {:.2f})".format(fdr_threshold))
    ax.set_xticks(x)
    ax.set_xticklabels(results_df.index, rotation=90)
    ax.legend()
    
    plt.tight_layout()
    bar_file = os.path.join(output_dir, f"{output_prefix}_bar.png")
    fig.savefig(bar_file, dpi=300)
    plt.close(fig)
    
    output_files["bar"] = bar_file
    logging.info(f"Saved essential gene enrichment bar chart to {bar_file}")
    
    # Save results as CSV
    csv_file = os.path.join(output_dir, f"{output_prefix}_summary.csv")
    results_df.to_csv(csv_file)
    output_files["summary"] = csv_file
    logging.info(f"Saved essential gene enrichment summary to {csv_file}")
    
    return output_files 