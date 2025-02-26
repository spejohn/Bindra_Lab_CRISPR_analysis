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
    run_drugz_analysis,
    find_design_matrix,
    process_mle_analysis
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
    skip_qc: bool = False,
    skip_mle: bool = False,
    use_docker: bool = True
) -> Dict[str, Any]:
    """
    Run the entire analysis pipeline on a directory of fastq files or count tables.
    
    Args:
        input_dir: Directory containing FASTQ files or count tables
        output_dir: Directory for output files
        library_file: Path to the library file
        experiment_name: Name for the experiment
        contrasts_file: Path to the contrasts file (for differential analysis)
        norm_method: Normalization method for MAGeCK
        fdr_threshold: FDR threshold for significant genes
        sample_sheet: Sample sheet for mapping sample names to conditions
        essential_genes: Dictionary of essential genes for QC
        overwrite: Whether to overwrite existing output files
        skip_drugz: Skip DrugZ analysis
        skip_qc: Skip quality control checks
        skip_mle: Skip MAGeCK MLE analysis
        use_docker: Use Docker containers for analysis tools when available
        
    Returns:
        Dictionary of results
    """
    # Check Docker availability if requested
    if use_docker:
        from analysis_pipeline.core.utils import check_docker_available
        docker_available = check_docker_available()
        if not docker_available:
            logging.warning("Docker is not available or required images are missing. Falling back to local installations.")
            use_docker = False
    
    # Create experiment directories
    dirs = create_experiment_dirs(output_dir, experiment_name)
    
    # Setup logging
    log_file = setup_logging(output_dir, experiment_name)
    
    # Log system information
    log_system_info()
    
    logging.info(f"Starting pipeline for {experiment_name}")
    logging.info(f"Input directory: {input_dir}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Library file: {library_file}")
    
    if contrasts_file:
        logging.info(f"Contrasts file: {contrasts_file}")
    
    # Initialize results dictionary
    results = {
        "experiment_name": experiment_name,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "library_file": library_file,
        "dirs": dirs,
        "log_file": log_file,
        "sample_sheet": None,
        "count_files": {},
        "rra_results": {},
        "mle_results": {},
        "drugz_results": {},
        "qc_results": {}
    }
    
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
    
    # Initialize a DataFrame to track read count files
    read_count_df = pd.DataFrame(columns=["file", "processed"])
    
    # Process read count files if present
    if input_files["read_counts"]:
        logging.info(f"Found {len(input_files['read_counts'])} read count files")
        
        for rc_file, info in input_files["read_counts"].items():
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
        
        for fastq_dir, info in input_files["fastq_dirs"].items():
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
            
            mageck_output = process_contrasts(
                contrasts_file=contrasts_file,
                count_table=count_file,
                output_dir=dirs["rra"],
                norm_method=norm_method,
                overwrite=overwrite
            )
            
            if "error" in mageck_output:
                logging.error(f"Error in MAGeCK analysis for {name}: {mageck_output['error']}")
            else:
                results["rra_results"][name] = mageck_output
                logging.info(f"Completed MAGeCK analysis for {name}")
            
            # Check if this experiment has a design matrix for MLE analysis
            contrast_dir = os.path.dirname(contrasts_file)
            design_matrix = None
            
            # First check if we already have a design matrix from input files
            fastq_dir_info = next((info for key, info in input_files["fastq_dirs"].items() if os.path.basename(key) == name), None)
            if fastq_dir_info and fastq_dir_info.get("design_matrix"):
                design_matrix = fastq_dir_info["design_matrix"]
            else:
                # Try to find a design matrix in the same directory as the contrasts file
                design_matrix = find_design_matrix(contrast_dir)
            
            # Run MLE analysis if design matrix is available and not skipped
            if design_matrix and not skip_mle:
                logging.info(f"Design matrix found for {name}, running MAGeCK MLE analysis")
                
                mle_output = process_mle_analysis(
                    count_table=count_file,
                    design_matrix=design_matrix,
                    output_dir=dirs["mle"],
                    experiment_name=f"{experiment_name}_{name}",
                    norm_method=norm_method,
                    overwrite=overwrite
                )
                
                if isinstance(mle_output, dict) and "error" in mle_output:
                    logging.error(f"Error in MAGeCK MLE analysis for {name}: {mle_output['error']}")
                else:
                    results["mle_results"][name] = mle_output
                    logging.info(f"Completed MAGeCK MLE analysis for {name}")
            elif not skip_mle:
                logging.info(f"No design matrix found for {name}, skipping MAGeCK MLE analysis")
    
    # Run DrugZ analysis if not skipped
    if not skip_drugz and contrasts_file and os.path.exists(contrasts_file):
        logging.info("Running DrugZ analysis...")
        
        for name, count_file in results["count_files"].items():
            drugz_results = run_drugz_analysis(
                count_table=count_file,
                contrasts_file=contrasts_file,
                output_dir=dirs["drugz"],
                overwrite=overwrite,
                use_docker=use_docker
            )
            
            if drugz_results and "error" not in drugz_results:
                results["drugz_results"][name] = drugz_results
                logging.info(f"Completed DrugZ analysis for {name}")
            else:
                error_msg = drugz_results.get("error", {}).get("message", "Unknown error")
                logging.error(f"Error in DrugZ analysis for {name}: {error_msg}")
    
    # Run QC plots if not skipped
    if not skip_qc:
        logging.info("Generating QC plots...")
        
        for name, count_file in results["count_files"].items():
            try:
                qc_plots = {}
                
                # Plot count distribution
                count_dist_plot = plot_count_distribution(
                    count_file=count_file,
                    output_dir=dirs["qc"],
                    output_prefix=f"{experiment_name}_{name}"
                )
                if count_dist_plot:
                    qc_plots["count_distribution"] = count_dist_plot
                
                # Plot guide correlations
                guide_corr_plot = plot_guide_correlations(
                    count_file=count_file,
                    output_dir=dirs["qc"],
                    output_prefix=f"{experiment_name}_{name}"
                )
                if guide_corr_plot:
                    qc_plots["guide_correlations"] = guide_corr_plot
                
                # Plot essential gene enrichment if available
                if essential_genes:
                    for gene_set_name, gene_list in essential_genes.items():
                        gene_enrichment_plot = plot_essential_gene_enrichment(
                            count_file=count_file,
                            gene_list=gene_list,
                            output_dir=dirs["qc"],
                            output_prefix=f"{experiment_name}_{name}_{gene_set_name}"
                        )
                        if gene_enrichment_plot:
                            qc_plots[f"essential_gene_enrichment_{gene_set_name}"] = gene_enrichment_plot
                
                results["qc_results"][name] = qc_plots
                
            except Exception as e:
                logging.error(f"Error generating QC plots for {name}: {e}")
    
    logging.info("Pipeline completed successfully")
    return results


def main():
    """Command-line entry point for the CRISPR screening analysis pipeline."""
    parser = argparse.ArgumentParser(description="CRISPR Screening Analysis Pipeline")
    
    # Required argument - just the input directory
    parser.add_argument("input_dir", help="Directory containing CRISPR screen data (FASTQ or count tables)")
    
    # Optional arguments
    parser.add_argument("-o", "--output-dir", help="Directory for output files (defaults to input_dir/results)")
    parser.add_argument("-l", "--library-file", help="Path to the library file (defaults to input_dir/library.csv)")
    parser.add_argument("-e", "--experiment-name", default="experiment", help="Name of the experiment")
    parser.add_argument("-c", "--contrasts-file", help="Path to the contrasts file (defaults to input_dir/contrasts.csv)")
    parser.add_argument("--norm-method", default=DEFAULT_NORM_METHOD, help="Normalization method")
    parser.add_argument("--fdr-threshold", type=float, default=DEFAULT_FDR_THRESHOLD, help="FDR threshold for significance")
    parser.add_argument("-s", "--sample-sheet", help="Path to the sample sheet (will be auto-generated if not provided)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip quality control plotting")
    parser.add_argument("--skip-mle", action="store_true", help="Skip MAGeCK MLE analysis")
    parser.add_argument("--use-docker", action="store_true", help="Use Docker containers for analysis tools when available")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Set default paths if not provided
    input_dir = os.path.abspath(args.input_dir)
    output_dir = args.output_dir or os.path.join(input_dir, "results")
    library_file = args.library_file or os.path.join(input_dir, "library.csv")
    contrasts_file = args.contrasts_file or os.path.join(input_dir, "contrasts.csv")
    
    # Check if required files exist
    if not os.path.exists(library_file):
        print(f"Error: Library file not found at {library_file}")
        print("Please provide a valid library file path with --library-file")
        return 1
    
    if args.contrasts_file is None and not os.path.exists(contrasts_file):
        print(f"Warning: No contrasts file found at {contrasts_file}")
        print("Differential analysis will be skipped.")
        contrasts_file = None
    
    # Set logging level
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(log_level)
    
    # Run the pipeline
    results = run_pipeline(
        input_dir=input_dir,
        output_dir=output_dir,
        library_file=library_file,
        experiment_name=args.experiment_name,
        contrasts_file=contrasts_file,
        norm_method=args.norm_method,
        fdr_threshold=args.fdr_threshold,
        sample_sheet=args.sample_sheet,
        overwrite=args.overwrite,
        skip_drugz=args.skip_drugz,
        skip_qc=args.skip_qc,
        skip_mle=args.skip_mle,
        use_docker=args.use_docker
    )
    
    if "error" in results:
        print(f"Error: {results['error']}")
        return 1
    
    print("Pipeline completed successfully")
    print(f"Log file: {results['log_file']}")
    print(f"Results directory: {output_dir}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 