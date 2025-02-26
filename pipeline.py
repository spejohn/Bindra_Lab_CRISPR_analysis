"""
CRISPR Screening Analysis Pipeline

This module provides the main entry point for the CRISPR screening analysis pipeline.
It ties together modules for sample processing, analysis, and quality control.
"""

import os
import sys
import argparse
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any

# Import core modules
from analysis_pipeline.core.config import (
    DEFAULT_NORM_METHOD,
    DEFAULT_ADJUST_METHOD,
    DEFAULT_FDR_THRESHOLD,
    DEFAULT_ESSENTIAL_GENES
)
from analysis_pipeline.core.logging_setup import setup_logging, log_system_info, log_input_parameters
from analysis_pipeline.core.file_handling import (
    ensure_output_dir,
    reverse_dir,
    find_input_files,
    create_experiment_dirs
)
from analysis_pipeline.core.validation import (
    validate_library_file,
    validate_input_tables,
    check_existing_files
)

# Import docker utilities
from analysis_pipeline.docker.docker_utils import verify_docker, test_docker_container

# Import analysis modules
from analysis_pipeline.analysis.sample_processing import (
    process_all_samples,
    merge_count_files,
    generate_sample_sheet
)
from analysis_pipeline.analysis.mageck_analysis import (
    process_contrasts,
    run_drugz_analysis
)

# Import QC modules
from analysis_pipeline.qc.plot_quality import (
    plot_count_distribution,
    plot_guide_correlations,
    plot_mageck_results,
    plot_essential_gene_enrichment
)


def run_pipeline(
    input_dir: str,
    output_dir: str,
    library_file: str,
    experiment_name: str = "experiment",
    contrasts_file: Optional[str] = None,
    norm_method: str = DEFAULT_NORM_METHOD,
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
    sample_sheet: Optional[str] = None,
    essential_genes: Optional[Dict[str, List[str]]] = None,
    overwrite: bool = False,
    skip_drugz: bool = False,
    skip_qc: bool = False
) -> Dict[str, Any]:
    """
    Run the CRISPR screening analysis pipeline.
    
    Args:
        input_dir: Directory containing input files (FASTQ or count tables)
        output_dir: Directory for output files
        library_file: Path to the library file
        experiment_name: Name of the experiment
        contrasts_file: Path to the contrasts file
        norm_method: Normalization method
        fdr_threshold: FDR threshold for significance
        sample_sheet: Path to the sample sheet
        essential_genes: Dictionary of essential gene lists
        overwrite: Whether to overwrite existing output files
        skip_drugz: Whether to skip DrugZ analysis
        skip_qc: Whether to skip quality control plotting
        
    Returns:
        Dictionary of results and output files
    """
    # Setup logging
    log_file = setup_logging(output_dir, experiment_name)
    
    # Log system information
    log_system_info()
    
    # Log input parameters
    log_input_parameters({
        "input_dir": input_dir,
        "output_dir": output_dir,
        "library_file": library_file,
        "experiment_name": experiment_name,
        "contrasts_file": contrasts_file,
        "norm_method": norm_method,
        "fdr_threshold": fdr_threshold,
        "sample_sheet": sample_sheet,
        "overwrite": overwrite,
        "skip_drugz": skip_drugz,
        "skip_qc": skip_qc
    })
    
    # Create experiment directories
    dirs = create_experiment_dirs(output_dir, experiment_name)
    
    # Verify Docker
    logging.info("Verifying Docker installation...")
    if not verify_docker():
        logging.error("Docker verification failed. Please check your Docker installation.")
        return {"error": "Docker verification failed"}
    
    # Test Docker container
    logging.info("Testing Docker container...")
    if not test_docker_container():
        logging.error("Docker container test failed. Please check your Docker installation and image.")
        return {"error": "Docker container test failed"}
    
    # Validate library file
    logging.info(f"Validating library file: {library_file}")
    valid_library = validate_library_file(library_file, dirs["library"])
    
    if not valid_library:
        logging.error("Library file validation failed. Please check your library file.")
        return {"error": "Library file validation failed"}
    
    # Check existing files
    logging.info("Checking for existing analyzed files...")
    existing_files = check_existing_files(input_dir, output_dir)
    logging.info(f"Found {len(existing_files)} previously analyzed files")
    
    # Find input files
    logging.info(f"Finding input files in {input_dir}")
    input_files = find_input_files(input_dir)
    
    # Initialize results
    results = {
        "experiment_name": experiment_name,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "library_file": library_file,
        "log_file": log_file,
        "directories": dirs,
        "count_files": {},
        "mageck_results": {},
        "drugz_results": {},
        "qc_plots": {}
    }
    
    # Initialize a DataFrame to track read count files
    read_count_df = pd.DataFrame(columns=["file", "processed"])
    
    # Process read count files if present
    if input_files["count_tables"]:
        logging.info(f"Found {len(input_files['count_tables'])} read count files")
        
        for rc_file in input_files["count_tables"]:
            if rc_file in existing_files and not overwrite:
                logging.info(f"Skipping previously processed read count file: {rc_file}")
                read_count_df = read_count_df.append({"file": rc_file, "processed": True}, ignore_index=True)
                continue
                
            try:
                # Process the read count file
                logging.info(f"Processing read count file: {rc_file}")
                # Add processing logic here
                
                read_count_df = read_count_df.append({"file": rc_file, "processed": True}, ignore_index=True)
            except Exception as e:
                logging.error(f"Error processing read count file {rc_file}: {e}")
                read_count_df = read_count_df.append({"file": rc_file, "processed": False, "error": str(e)}, ignore_index=True)
                
    # Initialize a DataFrame to track FASTQ files
    fastq_df = pd.DataFrame(columns=["file", "processed"])
    
    # Process FASTQ directories if present
    if input_files["fastq_dirs"]:
        logging.info(f"Found {len(input_files['fastq_dirs'])} FASTQ directories")
        
        for fastq_dir in input_files["fastq_dirs"]:
            # Create output directory for this FASTQ directory
            fastq_output_dir = reverse_dir(fastq_dir, input_dir, dirs["counts"])
            ensure_output_dir(fastq_output_dir)
            
            # Read or generate sample sheet
            sample_sheet_df = None
            
            if sample_sheet and os.path.exists(sample_sheet):
                try:
                    sample_sheet_df = pd.read_csv(sample_sheet)
                    logging.info(f"Read sample sheet from {sample_sheet} with {len(sample_sheet_df)} samples")
                except Exception as e:
                    logging.error(f"Error reading sample sheet {sample_sheet}: {e}")
            
            if sample_sheet_df is None:
                logging.info(f"Generating sample sheet for FASTQ directory: {fastq_dir}")
                sample_sheet_df = generate_sample_sheet(
                    fastq_dir=fastq_dir,
                    output_dir=dirs["samplesheets"],
                    experiment_name=f"{experiment_name}_{os.path.basename(fastq_dir)}"
                )
                
                if sample_sheet_df is None:
                    logging.error(f"No FASTQ files found in {fastq_dir}")
                    continue
            
            # Process all samples in the FASTQ directory
            count_files = process_all_samples(
                fastq_dir=fastq_dir,
                library_file=library_file,
                output_dir=fastq_output_dir,
                sample_sheet=sample_sheet_df,
                overwrite=overwrite
            )
            
            if not count_files:
                logging.error(f"No samples processed in {fastq_dir}")
                continue
            
            # Merge count files
            merged_file = merge_count_files(
                count_files=count_files,
                output_dir=dirs["counts"],
                output_name=f"{experiment_name}_{os.path.basename(fastq_dir)}"
            )
            
            if merged_file:
                results["count_files"][os.path.basename(fastq_dir)] = merged_file
                logging.info(f"Merged count file: {merged_file}")
            
    # Run MAGeCK analysis if contrasts file is provided
    if contrasts_file and os.path.exists(contrasts_file):
        logging.info(f"Running MAGeCK analysis with contrasts from {contrasts_file}")
        
        # Validate contrasts file
        validate_input_tables(contrasts_file)
        
        # For each merged count file, run MAGeCK analysis
        for name, count_file in results["count_files"].items():
            logging.info(f"Running MAGeCK analysis for count file: {count_file}")
            
            mageck_output_dir = os.path.join(dirs["mageck"], name)
            ensure_output_dir(mageck_output_dir)
            
            # Run MAGeCK analysis
            mageck_results = process_contrasts(
                contrasts_file=contrasts_file,
                count_table=count_file,
                output_dir=mageck_output_dir,
                norm_method=norm_method,
                overwrite=overwrite
            )
            
            results["mageck_results"][name] = mageck_results
            
            # Run DrugZ analysis if not skipped
            if not skip_drugz:
                logging.info(f"Running DrugZ analysis for count file: {count_file}")
                
                drugz_output_dir = os.path.join(dirs["drugz"], name)
                ensure_output_dir(drugz_output_dir)
                
                # Run DrugZ analysis
                drugz_results = run_drugz_analysis(
                    count_table=count_file,
                    contrasts_file=contrasts_file,
                    output_dir=drugz_output_dir,
                    overwrite=overwrite
                )
                
                results["drugz_results"][name] = drugz_results
    
    # Run quality control plots if not skipped
    if not skip_qc:
        logging.info("Generating quality control plots")
        
        qc_plots = {}
        
        # For each merged count file, generate count distribution plots
        for name, count_file in results["count_files"].items():
            logging.info(f"Generating count distribution plots for: {count_file}")
            
            count_plots_dir = os.path.join(dirs["qc"], f"{name}_counts")
            ensure_output_dir(count_plots_dir)
            
            # Plot count distribution
            count_plots = plot_count_distribution(
                count_table=count_file,
                output_dir=count_plots_dir,
                output_prefix=f"{name}_counts"
            )
            
            qc_plots[f"{name}_counts"] = count_plots
            
            # Plot guide correlations
            corr_plots_dir = os.path.join(dirs["qc"], f"{name}_correlations")
            ensure_output_dir(corr_plots_dir)
            
            corr_plots = plot_guide_correlations(
                count_table=count_file,
                output_dir=corr_plots_dir,
                output_prefix=f"{name}_correlations"
            )
            
            qc_plots[f"{name}_correlations"] = corr_plots
        
        # For each MAGeCK result, generate volcano plots
        for name, mageck_result in results["mageck_results"].items():
            if "error" in mageck_result:
                continue
                
            mageck_plots_dir = os.path.join(dirs["qc"], f"{name}_mageck")
            ensure_output_dir(mageck_plots_dir)
            
            # For each contrast, plot MAGeCK results
            for contrast, files in mageck_result.items():
                if "error" in files:
                    continue
                    
                if "gene_summary" not in files:
                    continue
                
                logging.info(f"Generating MAGeCK plots for contrast: {contrast}")
                
                mageck_plots = plot_mageck_results(
                    gene_summary_file=files["gene_summary"],
                    output_dir=mageck_plots_dir,
                    fdr_threshold=fdr_threshold,
                    essential_genes=essential_genes,
                    output_prefix=f"{name}_{contrast}"
                )
                
                qc_plots[f"{name}_{contrast}_mageck"] = mageck_plots
            
            # If there are multiple contrasts, plot essential gene enrichment
            gene_summary_files = {}
            for contrast, files in mageck_result.items():
                if "error" in files or "gene_summary" not in files:
                    continue
                gene_summary_files[contrast] = files["gene_summary"]
            
            if len(gene_summary_files) > 1:
                logging.info(f"Generating essential gene enrichment plot for {name}")
                
                enrichment_plots = plot_essential_gene_enrichment(
                    gene_summary_files=gene_summary_files,
                    output_dir=mageck_plots_dir,
                    essential_genes=essential_genes,
                    fdr_threshold=fdr_threshold,
                    output_prefix=f"{name}_enrichment"
                )
                
                qc_plots[f"{name}_enrichment"] = enrichment_plots
        
        results["qc_plots"] = qc_plots
    
    logging.info("Pipeline completed successfully")
    return results


def main():
    """Command-line entry point for the CRISPR screening analysis pipeline."""
    parser = argparse.ArgumentParser(description="CRISPR Screening Analysis Pipeline")
    
    # Required arguments
    parser.add_argument("--input-dir", required=True, help="Directory containing input files (FASTQ or count tables)")
    parser.add_argument("--output-dir", required=True, help="Directory for output files")
    parser.add_argument("--library-file", required=True, help="Path to the library file")
    
    # Optional arguments
    parser.add_argument("--experiment-name", default="experiment", help="Name of the experiment")
    parser.add_argument("--contrasts-file", help="Path to the contrasts file")
    parser.add_argument("--norm-method", default=DEFAULT_NORM_METHOD, help="Normalization method")
    parser.add_argument("--fdr-threshold", type=float, default=DEFAULT_FDR_THRESHOLD, help="FDR threshold for significance")
    parser.add_argument("--sample-sheet", help="Path to the sample sheet")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip quality control plotting")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Run the pipeline
    results = run_pipeline(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        library_file=args.library_file,
        experiment_name=args.experiment_name,
        contrasts_file=args.contrasts_file,
        norm_method=args.norm_method,
        fdr_threshold=args.fdr_threshold,
        sample_sheet=args.sample_sheet,
        overwrite=args.overwrite,
        skip_drugz=args.skip_drugz,
        skip_qc=args.skip_qc
    )
    
    if "error" in results:
        print(f"Error: {results['error']}")
        sys.exit(1)
    
    print("Pipeline completed successfully")
    print(f"Log file: {results['log_file']}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 