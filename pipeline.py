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
import glob

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
    create_experiment_dirs,
    convert_results_to_csv
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


def process_experiment_by_type(
    input_dir: str,
    output_dir: str,
    library_file: str,
    experiment_name: str,
    contrasts_file: Optional[str] = None,
    contrasts: Optional[List[str]] = None,
    norm_method: str = DEFAULT_NORM_METHOD,
    sample_sheet: Optional[str] = None,
    overwrite: bool = False,
    skip_drugz: bool = False,
    skip_mle: bool = False,
    use_docker: bool = True,
    data_type: str = "fastq"  # Either "fastq" or "count"
) -> Dict[str, Any]:
    """
    Process an experiment with a specific data type (FASTQ or count).
    
    Args:
        input_dir: Directory containing FASTQ files or count tables
        output_dir: Directory for output files
        library_file: Path to the library file
        experiment_name: Name for the experiment (with suffix)
        contrasts_file: Path to the contrasts file
        contrasts: List of contrast names
        norm_method: Normalization method for MAGeCK
        sample_sheet: Sample sheet for mapping sample names to conditions
        overwrite: Whether to overwrite existing output files
        skip_drugz: Skip DrugZ analysis
        skip_mle: Skip MAGeCK MLE analysis
        use_docker: Use Docker containers for analysis tools when available
        data_type: Type of data to process (fastq or count)
        
    Returns:
        Dictionary of results
    """
    # Create base experiment directory
    base_exp_dir = os.path.join(output_dir, experiment_name)
    os.makedirs(base_exp_dir, exist_ok=True)
    
    # Setup a single log file for the entire experiment
    log_file = os.path.join(base_exp_dir, f"{experiment_name}.log")
    setup_logging(log_file=log_file)
    
    # Log system information
    log_system_info()
    
    logging.info(f"Starting pipeline for {experiment_name} with {data_type} data")
    logging.info(f"Input directory: {input_dir}")
    logging.info(f"Output directory: {base_exp_dir}")
    logging.info(f"Library file: {library_file}")
    
    if contrasts_file:
        logging.info(f"Contrasts file: {contrasts_file}")
    
    # Initialize results dictionary
    results = {
        "experiment_name": experiment_name,
        "input_dir": input_dir,
        "output_dir": base_exp_dir,
        "library_file": library_file,
        "log_file": log_file,
        "contrasts": {},
        "count_files": {}
    }
    
    # Generate sample sheet for the experiment
    if sample_sheet is None:
        sample_sheet_path = os.path.join(base_exp_dir, f"{experiment_name}_samples.txt")
        generate_sample_sheet(input_dir, base_exp_dir, experiment_name, experiment_name=experiment_name)
        results["sample_sheet"] = sample_sheet_path
    else:
        # Copy the user-provided sample sheet to the experiment directory
        sample_sheet_path = os.path.join(base_exp_dir, f"{experiment_name}_samples.txt")
        import shutil
        shutil.copy(sample_sheet, sample_sheet_path)
        results["sample_sheet"] = sample_sheet_path
    
    count_table = None
    
    if data_type == "fastq":
        # Process all FASTQ samples
        fastq_dirs = glob.glob(os.path.join(input_dir, '**/fastq'), recursive=True)
        logging.info(f"Processing {len(fastq_dirs)} FASTQ directories")
        
        count_files = process_all_samples(input_dir, library_file, base_exp_dir, sample_sheet=pd.read_csv(sample_sheet_path), overwrite=overwrite, experiment_name=experiment_name)
        
        # Merge count files into a single table
        count_table = merge_count_files(count_files, base_exp_dir, experiment_name)
        if count_table:
            results["count_table"] = count_table
            logging.info(f"Created merged count table: {count_table}")
        else:
            logging.error("Failed to create merged count table")
            return {"error": "Failed to create merged count table"}
    
    elif data_type == "count":
        # Look for existing count tables
        count_tables = glob.glob(os.path.join(input_dir, '**/*.count'), recursive=True)
        if count_tables:
            count_table = count_tables[0]
            # Copy the count table to the experiment directory
            dest_path = os.path.join(base_exp_dir, f"{experiment_name}.count")
            import shutil
            shutil.copy(count_table, dest_path)
            count_table = dest_path
            results["count_table"] = count_table
            logging.info(f"Using existing count table: {count_table}")
        else:
            logging.error("No count tables found")
            return {"error": "No count tables found"}
    
    # Run analysis for each contrast
    if contrasts and count_table:
        for contrast in contrasts:
            logging.info(f"Processing contrast: {contrast}")
            
            # Create a directory for this contrast
            contrast_dir = os.path.join(base_exp_dir, contrast)
            os.makedirs(contrast_dir, exist_ok=True)
            
            # Initialize results for this contrast
            results["contrasts"][contrast] = {
                "dir": contrast_dir
            }
            
            # Run RRA analysis
            if contrasts_file:
                # Run RRA for this contrast
                rra_results = process_contrasts(
                    contrasts_file=contrasts_file,
                    count_table=count_table,
                    output_dir=contrast_dir,
                    norm_method=norm_method,
                    overwrite=overwrite,
                    target_contrast=contrast
                )
                
                # Convert RRA results to CSV
                if rra_results and contrast in rra_results and "gene_summary" in rra_results[contrast]:
                    rra_csv = convert_results_to_csv(rra_results[contrast]["gene_summary"], "RRA")
                    rra_results[contrast]["csv"] = rra_csv
                
                results["contrasts"][contrast]["rra_results"] = rra_results
            
            # Run MLE analysis if design matrix is available and not skipped
            if not skip_mle:
                design_matrix = find_design_matrix(input_dir)
                
                if design_matrix:
                    logging.info(f"Found design matrix: {design_matrix}")
                    
                    # Run MLE analysis
                    mle_results = process_mle_analysis(
                        count_table=count_table,
                        design_matrix=design_matrix,
                        output_dir=contrast_dir,
                        contrast_name=contrast,
                        norm_method=norm_method,
                        overwrite=overwrite
                    )
                    
                    # Convert MLE results to CSV
                    if mle_results and "gene_summary" in mle_results:
                        mle_csv = convert_results_to_csv(mle_results["gene_summary"], "MLE")
                        mle_results["csv"] = mle_csv
                    
                    results["contrasts"][contrast]["mle_results"] = mle_results
                else:
                    logging.warning(f"No design matrix found, skipping MLE analysis for {contrast}")
            
            # Run DrugZ analysis if requested
            if not skip_drugz:
                logging.info(f"Running DrugZ analysis for {contrast}")
                
                drugz_results = run_drugz_analysis(
                    count_table=count_table,
                    contrasts_file=contrasts_file,
                    output_dir=contrast_dir,
                    overwrite=overwrite,
                    use_docker=use_docker,
                    target_contrast=contrast
                )
                
                # Convert DrugZ results to CSV
                if drugz_results and contrast in drugz_results and "drugz_scores" in drugz_results[contrast]:
                    drugz_csv = convert_results_to_csv(drugz_results[contrast]["drugz_scores"], "DrugZ")
                    drugz_results[contrast]["csv"] = drugz_csv
                
                results["contrasts"][contrast]["drugz_results"] = drugz_results
    
    return results


def run_pipeline(
    input_dir: str,
    output_dir: str,
    library_file: Optional[str] = None,
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
        input_dir: Directory containing experiment directories
        output_dir: Directory for output files
        library_file: Path to the library file (if not provided, will look in experiment directory)
        experiment_name: Name of the experiment (corresponds to subdirectory in input_dir)
        contrasts_file: Path to the contrasts file (if not provided, will look in experiment directory)
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
    
    # Look for input files using experiment_name as subdirectory
    from analysis_pipeline.core.file_handling import find_input_files
    input_files = find_input_files(input_dir, experiment_name)
    
    # If library_file is not provided, check if one was found in the experiment directory
    if library_file is None and input_files.get('library_file'):
        library_file = input_files['library_file']
        logging.info(f"Using library file from experiment directory: {library_file}")
    elif library_file is None:
        # Use default library location
        library_file = os.path.join(input_dir, "library.csv")
        
    # If contrasts_file is not provided, check if one was found in the experiment directory
    if contrasts_file is None and input_files.get('contrast_file'):
        contrasts_file = input_files['contrast_file']
        logging.info(f"Using contrast file from experiment directory: {contrasts_file}")
    elif contrasts_file is None:
        # Use default contrasts location
        contrasts_file = os.path.join(input_dir, "contrasts.csv")
    
    # Read contrasts file to get contrast names
    contrasts = []
    if contrasts_file and os.path.exists(contrasts_file):
        try:
            contrasts_df = pd.read_csv(contrasts_file)
            if "contrast" in contrasts_df.columns:
                contrasts = contrasts_df["contrast"].tolist()
                logging.info(f"Found {len(contrasts)} contrasts: {', '.join(contrasts)}")
            else:
                logging.error("Contrasts file does not contain a 'contrast' column")
        except Exception as e:
            logging.error(f"Error reading contrasts file: {e}")
    
    # Initialize results dictionary
    results = {
        "experiment_name": experiment_name,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "library_file": library_file,
        "contrasts": {},
        "count_files": {}
    }
    
    # Determine input data types - look directly in experiment subdirectory or data dir
    # Check in experiment directory first
    fastq_pattern = os.path.join(input_dir, experiment_name, "**", "*.fastq*")
    fastq_files = glob.glob(fastq_pattern, recursive=True)
    
    # If no fastq files found in experiment dir, check in data dir
    if not fastq_files:
        fastq_pattern = os.path.join(input_dir, experiment_name, "data", "**", "*.fastq*")
        fastq_files = glob.glob(fastq_pattern, recursive=True)
    
    # Also look for fastq in fastq subdirectory
    if not fastq_files:
        fastq_pattern = os.path.join(input_dir, experiment_name, "fastq", "**", "*.fastq*")
        fastq_files = glob.glob(fastq_pattern, recursive=True)
        
    # Also try fq extension
    if not fastq_files:
        fastq_pattern = os.path.join(input_dir, experiment_name, "**", "*.fq*")
        fastq_files = glob.glob(fastq_pattern, recursive=True)
        
    # Check count files
    count_pattern = os.path.join(input_dir, experiment_name, "**", "*.count")
    count_files = glob.glob(count_pattern, recursive=True)
    
    # If no count files found, check in data dir
    if not count_files:
        count_pattern = os.path.join(input_dir, experiment_name, "data", "**", "*.count")
        count_files = glob.glob(count_pattern, recursive=True)
    
    # Also look for count files in counts subdirectory
    if not count_files:
        count_pattern = os.path.join(input_dir, experiment_name, "counts", "**", "*.count")
        count_files = glob.glob(count_pattern, recursive=True)
        
    has_fastq = len(fastq_files) > 0
    has_count_files = len(count_files) > 0
    
    logging.info(f"Found {len(fastq_files)} FASTQ files and {len(count_files)} count files")
    
    # If both types exist, create separate experiment directories
    if has_fastq and has_count_files:
        logging.info(f"Found both FASTQ files and pre-processed count files. Creating separate analysis directories.")
        
        # Process FASTQ files
        if has_fastq:
            fastq_exp_name = f"{experiment_name}_FASTQ"
            fastq_results = process_experiment_by_type(
                input_dir=input_dir,
                output_dir=output_dir,
                library_file=library_file,
                experiment_name=fastq_exp_name,
                contrasts_file=contrasts_file,
                contrasts=contrasts,
                norm_method=norm_method,
                sample_sheet=sample_sheet,
                overwrite=overwrite,
                skip_drugz=skip_drugz,
                skip_mle=skip_mle,
                use_docker=use_docker,
                data_type="fastq"
            )
            results["fastq_analysis"] = fastq_results
        
        # Process count files
        if has_count_files:
            count_exp_name = f"{experiment_name}_RC"
            count_results = process_experiment_by_type(
                input_dir=input_dir,
                output_dir=output_dir,
                library_file=library_file,
                experiment_name=count_exp_name,
                contrasts_file=contrasts_file,
                contrasts=contrasts,
                norm_method=norm_method,
                sample_sheet=sample_sheet,
                overwrite=overwrite,
                skip_drugz=skip_drugz,
                skip_mle=skip_mle,
                use_docker=use_docker,
                data_type="count"
            )
            results["count_analysis"] = count_results
        
        return results
    else:
        # Only one type exists, proceed with standard processing
        # Create base experiment directory
        base_exp_dir = os.path.join(output_dir, experiment_name)
        os.makedirs(base_exp_dir, exist_ok=True)
        
        # Setup a single log file for the entire experiment
        log_file = os.path.join(base_exp_dir, f"{experiment_name}.log")
        setup_logging(log_file=log_file)
        
        # Log system information
        log_system_info()
        
        logging.info(f"Starting pipeline for {experiment_name}")
        logging.info(f"Input directory: {input_dir}")
        logging.info(f"Output directory: {output_dir}")
        logging.info(f"Library file: {library_file}")
        
        if contrasts_file:
            logging.info(f"Contrasts file: {contrasts_file}")
        
        results["log_file"] = log_file
        
        # Generate sample sheet for the experiment
        if sample_sheet is None:
            sample_sheet_path = os.path.join(base_exp_dir, f"{experiment_name}_samples.txt")
            generate_sample_sheet(input_dir, base_exp_dir, experiment_name, experiment_name=experiment_name)
            results["sample_sheet"] = sample_sheet_path
        else:
            # Copy the user-provided sample sheet to the experiment directory
            sample_sheet_path = os.path.join(base_exp_dir, f"{experiment_name}_samples.txt")
            import shutil
            shutil.copy(sample_sheet, sample_sheet_path)
            results["sample_sheet"] = sample_sheet_path
        
        count_table = None
        
        if has_fastq:
            logging.info(f"Found {len(fastq_files)} fastq files")
            # Process all samples
            count_files = process_all_samples(
                input_dir, 
                library_file, 
                base_exp_dir, 
                sample_sheet=pd.read_csv(sample_sheet_path), 
                overwrite=overwrite,
                experiment_name=experiment_name
            )
            
            # Merge count files into a single table - save directly in the experiment directory
            count_table = merge_count_files(count_files, base_exp_dir, experiment_name)
            if count_table:
                results["count_table"] = count_table
                logging.info(f"Created merged count table: {count_table}")
            else:
                logging.error("Failed to create merged count table")
                return {"error": "Failed to create merged count table"}
        else:
            # Look for existing count tables
            count_tables = glob.glob(os.path.join(input_dir, '**/*.count'), recursive=True)
            if count_tables:
                count_table = count_tables[0]
                # Copy the count table to the experiment directory
                dest_path = os.path.join(base_exp_dir, f"{experiment_name}.count")
                import shutil
                shutil.copy(count_table, dest_path)
                count_table = dest_path
                results["count_table"] = count_table
                logging.info(f"Using existing count table: {count_table}")
            else:
                logging.error("No fastq files or count tables found")
                return {"error": "No fastq files or count tables found"}
        
        # Run analysis for each contrast
        if contrasts and count_table:
            for contrast in contrasts:
                logging.info(f"Processing contrast: {contrast}")
                
                # Create a directory for this contrast
                contrast_dir = os.path.join(base_exp_dir, contrast)
                os.makedirs(contrast_dir, exist_ok=True)
                
                # Initialize results for this contrast
                results["contrasts"][contrast] = {
                    "dir": contrast_dir
                }
                
                # Run RRA analysis
                if contrasts_file:
                    # Run RRA for this contrast
                    rra_results = process_contrasts(
                        contrasts_file=contrasts_file,
                        count_table=count_table,  # Use the main count table
                        output_dir=contrast_dir,   # Save directly in contrast directory
                        norm_method=norm_method,
                        overwrite=overwrite,
                        target_contrast=contrast
                    )
                    
                    # Convert RRA results to CSV
                    if rra_results and contrast in rra_results and "gene_summary" in rra_results[contrast]:
                        rra_csv = convert_results_to_csv(rra_results[contrast]["gene_summary"], "RRA")
                        rra_results[contrast]["csv"] = rra_csv
                    
                    results["contrasts"][contrast]["rra_results"] = rra_results
                
                # Run MLE analysis if design matrix is available and not skipped
                if not skip_mle:
                    design_matrix = find_design_matrix(input_dir)
                    
                    if design_matrix:
                        logging.info(f"Found design matrix: {design_matrix}")
                        
                        # Run MLE analysis
                        mle_results = process_mle_analysis(
                            count_table=count_table,  # Use the main count table
                            design_matrix=design_matrix,
                            output_dir=contrast_dir,   # Save directly in contrast directory
                            contrast_name=contrast,
                            norm_method=norm_method,
                            overwrite=overwrite
                        )
                        
                        # Convert MLE results to CSV
                        if mle_results and "gene_summary" in mle_results:
                            mle_csv = convert_results_to_csv(mle_results["gene_summary"], "MLE")
                            mle_results["csv"] = mle_csv
                        
                        results["contrasts"][contrast]["mle_results"] = mle_results
                    else:
                        logging.warning(f"No design matrix found, skipping MLE analysis for {contrast}")
                
                # Run DrugZ analysis if requested
                if not skip_drugz:
                    logging.info(f"Running DrugZ analysis for {contrast}")
                    
                    drugz_results = run_drugz_analysis(
                        count_table=count_table,  # Use the main count table
                        contrasts_file=contrasts_file,
                        output_dir=contrast_dir,   # Save directly in contrast directory
                        overwrite=overwrite,
                        use_docker=use_docker,
                        target_contrast=contrast
                    )
                    
                    # Convert DrugZ results to CSV
                    if drugz_results and contrast in drugz_results and "drugz_scores" in drugz_results[contrast]:
                        drugz_csv = convert_results_to_csv(drugz_results[contrast]["drugz_scores"], "DrugZ")
                        drugz_results[contrast]["csv"] = drugz_csv
                    
                    results["contrasts"][contrast]["drugz_results"] = drugz_results
        else:
            logging.warning("No contrasts or count table available, skipping differential analysis")
        
        logging.info(f"Analysis pipeline for {experiment_name} completed successfully")
        
        return results


def main():
    """Command-line entry point for the CRISPR screening analysis pipeline."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="CRISPR screening analysis pipeline")
    
    # Input/output options
    parser.add_argument("--input-dir", "-i", required=True, help="Directory containing input files or experiment directories")
    parser.add_argument("--output-dir", "-o", required=True, help="Directory for output files")
    parser.add_argument("--experiment-name", "-e", default="experiment", 
                        help="Name of the experiment; corresponds to a subdirectory in input-dir")
    parser.add_argument("--library-file", "-l", help="Path to library file (defaults to <input-dir>/<experiment-name>/library.csv)")
    parser.add_argument("--contrasts-file", "-c", help="Path to contrasts file (defaults to <input-dir>/<experiment-name>/contrasts.csv)")
    parser.add_argument("--sample-sheet", "-s", help="Path to sample sheet (will be generated if not provided)")
    
    # Analysis options
    parser.add_argument("--norm-method", default=DEFAULT_NORM_METHOD, 
                        choices=["median", "total", "control"], 
                        help="Normalization method for MAGeCK")
    parser.add_argument("--fdr-threshold", type=float, default=DEFAULT_FDR_THRESHOLD, 
                        help="FDR threshold for significant genes")
    
    # Flags
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip quality control checks")
    parser.add_argument("--skip-mle", action="store_true", help="Skip MAGeCK MLE analysis")
    parser.add_argument("--use-docker", action="store_true", help="Use Docker containers when available")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Create absolute paths
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Default library and contrasts files within the experiment directory
    library_file = args.library_file or os.path.join(input_dir, args.experiment_name, "library.csv")
    contrasts_file = args.contrasts_file or os.path.join(input_dir, args.experiment_name, "contrasts.csv")
    
    # Check if required files exist
    if not os.path.exists(library_file):
        # Also check legacy location
        legacy_library = os.path.join(input_dir, "library.csv")
        if os.path.exists(legacy_library):
            library_file = legacy_library
            print(f"Using legacy library file: {library_file}")
        else:
            print(f"Error: Library file not found at {library_file} or {legacy_library}")
            print("Please provide a valid library file path with --library-file")
            return 1
    
    if args.contrasts_file is None and not os.path.exists(contrasts_file):
        # Also check legacy location
        legacy_contrasts = os.path.join(input_dir, "contrasts.csv")
        if os.path.exists(legacy_contrasts):
            contrasts_file = legacy_contrasts
            print(f"Using legacy contrasts file: {contrasts_file}")
        else:
            print(f"Warning: No contrasts file found at {contrasts_file} or {legacy_contrasts}")
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