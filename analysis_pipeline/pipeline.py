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
import concurrent.futures
import time

# Import core modules
from analysis_pipeline.core.config import (
    DEFAULT_NORM_METHOD,
    DEFAULT_ADJUST_METHOD,
    DEFAULT_FDR_THRESHOLD,
    DEFAULT_ESSENTIAL_GENES,
    COUNT_CSV_PATTERN
)
from analysis_pipeline.core.logging_setup import setup_logging, log_system_info, log_input_parameters, ProgressReporter
from analysis_pipeline.core.file_handling import (
    ensure_output_dir,
    reverse_dir,
    find_input_files,
    create_experiment_dirs,
    convert_results_to_csv,
    make_count_table,
    copy_file,
    ensure_count_formats
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
    use_docker: bool = True,  # Always use Docker by default
    data_type: str = "fastq",  # Either "fastq" or "count"
    parallel: bool = True,
    max_workers: int = 4
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
        use_docker: Use Docker containers for analysis tools (always True, required for analysis)
        data_type: Type of data to process (fastq or count)
        parallel: Whether to use parallel processing for sample counting
        max_workers: Maximum number of parallel workers
        
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
        
        count_files = process_all_samples(
            input_dir, 
            library_file, 
            base_exp_dir, 
            sample_sheet=pd.read_csv(sample_sheet_path), 
            overwrite=overwrite, 
            experiment_name=experiment_name,
            parallel=parallel,
            max_workers=max_workers
        )
        
        # Merge count files into a single table
        count_table = merge_count_files(count_files, base_exp_dir, experiment_name, from_fastq=True)
        if count_table:
            results["count_table"] = count_table
            logging.info(f"Created merged count table from FASTQ files: {count_table}")
        else:
            logging.error("Failed to create merged count table")
            return {"error": "Failed to create merged count table"}
    
    elif data_type == "count":
        # Look for existing count tables or CSV files
        count_tables = glob.glob(os.path.join(input_dir, '**/*.count'), recursive=True)
        csv_tables = glob.glob(os.path.join(input_dir, '**/*.csv'), recursive=True)
        
        # Filter and prioritize CSV tables that look like count data
        count_related_csv = []
        for csv_file in csv_tables:
            base_name = os.path.basename(csv_file).lower()
            if 'count' in base_name or 'read_count' in base_name or '_rc' in base_name:
                count_related_csv.append(csv_file)
                
        # Log what we found
        if count_related_csv:
            logging.info(f"Found {len(count_related_csv)} CSV files that appear to be count data")
            
        if count_tables:
            count_table = count_tables[0]
            # Copy the count table to the experiment directory
            dest_path = os.path.join(base_exp_dir, f"{experiment_name}.count")
            import shutil
            shutil.copy(count_table, dest_path)
            count_table = dest_path
            results["count_table"] = count_table
            logging.info(f"Using existing count table: {count_table}")
        elif count_related_csv:
            # Use the first count-related CSV file
            csv_file = count_related_csv[0]
            from analysis_pipeline.core.file_handling import ensure_count_formats
            csv_path, tab_path = ensure_count_formats(csv_file, base_exp_dir)
            count_table = tab_path  # Use the tab-delimited file for analysis
            results["count_table"] = count_table
            logging.info(f"Created both formats - CSV: {csv_path}, Tab-delimited: {tab_path}")
        elif csv_tables:
            # If no count-specific CSV files found, use the first available CSV
            csv_file = csv_tables[0]
            from analysis_pipeline.core.file_handling import ensure_count_formats
            csv_path, tab_path = ensure_count_formats(csv_file, base_exp_dir)
            count_table = tab_path  # Use the tab-delimited file for analysis
            results["count_table"] = count_table
            logging.info(f"Created both formats - CSV: {csv_path}, Tab-delimited: {tab_path}")
        else:
            logging.error("No count tables or CSV files found")
            return {"error": "No count tables or CSV files found"}
    
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


def process_contrast_parallel(
    contrast: str,
    count_table: str,
    base_exp_dir: str,
    contrasts_file: Optional[str],
    norm_method: str,
    overwrite: bool,
    skip_drugz: bool,
    skip_mle: bool,
    use_docker: bool
) -> Dict[str, Any]:
    """
    Process a single contrast in parallel execution.
    
    Args:
        contrast: Name of the contrast to process
        count_table: Path to the count table
        base_exp_dir: Base experiment directory
        contrasts_file: Path to the contrasts file
        norm_method: Normalization method
        overwrite: Whether to overwrite existing files
        skip_drugz: Whether to skip DrugZ analysis
        skip_mle: Whether to skip MLE analysis
        use_docker: Whether to use Docker
        
    Returns:
        Dictionary with results for this contrast
    """
    try:
        logging.info(f"Processing contrast in parallel thread: {contrast}")
        
        # Create a directory for this contrast
        contrast_dir = os.path.join(base_exp_dir, contrast)
        os.makedirs(contrast_dir, exist_ok=True)
        
        # Initialize results for this contrast
        contrast_results = {
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
            
            contrast_results["rra_results"] = rra_results
        
        # Run MLE analysis if design matrix is available and not skipped
        if not skip_mle:
            design_matrix = find_design_matrix(os.path.dirname(count_table))
            
            if design_matrix:
                logging.info(f"Found design matrix for contrast {contrast}: {design_matrix}")
                
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
                
                contrast_results["mle_results"] = mle_results
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
            
            contrast_results["drugz_results"] = drugz_results
        
        return {contrast: contrast_results}
    
    except Exception as e:
        logging.error(f"Error processing contrast {contrast}: {str(e)}")
        return {contrast: {"error": str(e)}}


def run_analysis_contrasts_parallel(
    contrasts: List[str],
    count_table: str,
    base_exp_dir: str,
    contrasts_file: Optional[str],
    norm_method: str,
    overwrite: bool = False,
    skip_drugz: bool = False,
    skip_mle: bool = False,
    use_docker: bool = True,
    max_workers: int = 4
) -> Dict[str, Dict[str, Any]]:
    """
    Run analysis for multiple contrasts in parallel.
    
    Args:
        contrasts: List of contrast names
        count_table: Path to the count table
        base_exp_dir: Base experiment directory
        contrasts_file: Path to the contrasts file
        norm_method: Normalization method
        overwrite: Whether to overwrite existing files
        skip_drugz: Whether to skip DrugZ analysis
        skip_mle: Whether to skip MLE analysis
        use_docker: Whether to use Docker
        max_workers: Maximum number of parallel workers
        
    Returns:
        Dictionary with results for each contrast
    """
    if not contrasts:
        return {}
    
    # Determine number of workers (use min of contrasts or max_workers)
    n_workers = min(len(contrasts), max_workers)
    logging.info(f"Processing {len(contrasts)} contrasts using {n_workers} parallel workers")
    
    all_results = {}
    
    # Use ThreadPoolExecutor for parallel processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        # Submit all contrasts for processing
        future_to_contrast = {}
        for contrast in contrasts:
            future = executor.submit(
                process_contrast_parallel,
                contrast=contrast,
                count_table=count_table,
                base_exp_dir=base_exp_dir,
                contrasts_file=contrasts_file,
                norm_method=norm_method,
                overwrite=overwrite,
                skip_drugz=skip_drugz,
                skip_mle=skip_mle,
                use_docker=use_docker
            )
            future_to_contrast[future] = contrast
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_contrast):
            contrast = future_to_contrast[future]
            try:
                result = future.result()
                all_results.update(result)
                logging.info(f"Completed processing contrast: {contrast}")
            except Exception as e:
                logging.error(f"Error in parallel processing for contrast {contrast}: {str(e)}")
                all_results[contrast] = {"error": str(e)}
    
    return all_results


def run_pipeline(
    input_dir: str,
    output_dir: str,
    library_file: Optional[str] = None,
    experiment_name: Optional[str] = None,  # Changed to Optional with None default
    contrasts_file: Optional[str] = None,
    norm_method: str = DEFAULT_NORM_METHOD,
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
    sample_sheet: Optional[str] = None,
    essential_genes: Optional[Dict[str, List[str]]] = None,
    overwrite: bool = False,  # Default is False
    skip_drugz: bool = False,
    skip_qc: bool = False,
    skip_mle: bool = False,
    use_docker: bool = True,
    count_file: Optional[str] = None,
    parallel: bool = True,
    max_workers: int = 4,
    progress_reporter: Optional['ProgressReporter'] = None
) -> Dict[str, Any]:
    """
    Run the entire analysis pipeline on a directory of fastq files or count tables.
    
    Args:
        input_dir: Directory containing experiment directories
        output_dir: Directory for output files
        library_file: Path to the library file (if not provided, will look in experiment directory)
        experiment_name: Name of the experiment (if not provided, will use the directory name)
        contrasts_file: Path to the contrasts file (if not provided, will look in experiment directory)
        norm_method: Normalization method for MAGeCK
        fdr_threshold: FDR threshold for significant genes
        sample_sheet: Sample sheet for mapping sample names to conditions
        essential_genes: Dictionary of essential genes for QC
        overwrite: Whether to overwrite existing output files (defaults to False)
        skip_drugz: Skip DrugZ analysis
        skip_qc: Skip quality control checks
        skip_mle: Skip MAGeCK MLE analysis
        use_docker: Use Docker containers for analysis tools when available
        count_file: Path to a count file (.count or .csv) to use directly instead of looking in input-dir
        parallel: Whether to use parallel processing for sample counting and contrast analysis
        max_workers: Maximum number of parallel workers
        progress_reporter: Optional progress reporter instance
        
    Returns:
        Dictionary of results
    """
    # If experiment_name is not provided, derive it from the input directory
    if experiment_name is None:
        experiment_name = os.path.basename(os.path.normpath(input_dir))
        logging.info(f"No experiment name provided, using directory name: {experiment_name}")
    
    # Create a new progress reporter if none was provided
    if progress_reporter is None:
        # Import here to avoid circular imports
        from analysis_pipeline.core.logging_setup import ProgressReporter
        # Estimate steps based on available information
        total_steps = 10  # Default estimate
        progress_reporter = ProgressReporter(total_steps, experiment_name)
    
    # Check Docker availability if requested
    if use_docker:
        progress_reporter.update("Docker Check", "Checking Docker availability")
        from analysis_pipeline.core.utils import check_docker_available
        docker_available = check_docker_available()
        if not docker_available:
            logging.warning("Docker is not available or required images are missing. Falling back to local installations.")
            use_docker = False
    
    # Look for input files using experiment_name as subdirectory
    progress_reporter.update("Input Files", f"Scanning input directory: {input_dir}")
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
        progress_reporter.update("Contrasts", "Reading contrasts from file")
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
    
    # If a specific count file was provided, use it directly
    if count_file and os.path.exists(count_file):
        progress_reporter.update("Count File", f"Processing user-provided count file: {count_file}")
        logging.info(f"Using user-provided count file: {count_file}")
        # Copy the count file to the experiment directory
        base_exp_dir = os.path.join(output_dir, experiment_name)
        os.makedirs(base_exp_dir, exist_ok=True)
        
        # If it's a CSV file, convert it to tab-delimited format
        if count_file.endswith('.csv'):
            try:
                from analysis_pipeline.core.file_handling import ensure_count_formats
                progress.update("CSV Conversion", f"Converting CSV to tab-delimited: {count_file}")
                csv_path, tab_path = ensure_count_formats(count_file, output_dir)
                count_file = tab_path  # Use the tab-delimited file for analysis
                print(f"Created both formats - CSV: {csv_path}, Tab-delimited: {tab_path}")
            except Exception as e:
                print(f"Error converting count file formats: {str(e)}")
                return 1
        # If it's a tab-delimited file, ensure we also have a CSV version
        elif count_file.endswith('.count') or count_file.endswith('.counts') or count_file.endswith('.txt'):
            try:
                from analysis_pipeline.core.file_handling import ensure_count_formats
                progress.update("Format Conversion", f"Ensuring both formats exist for: {count_file}")
                csv_path, tab_path = ensure_count_formats(count_file, output_dir)
                count_file = tab_path  # Continue using the tab-delimited file for analysis
                print(f"Created both formats - CSV: {csv_path}, Tab-delimited: {tab_path}")
            except Exception as e:
                print(f"Error ensuring count file formats: {str(e)}")
                # Continue with the original file
        
        # Use our enhanced copy function with retry
        try:
            copy_file(count_file, os.path.join(base_exp_dir, f"{experiment_name}.count"))
            count_table = os.path.join(base_exp_dir, f"{experiment_name}.count")
            results["count_table"] = count_table
            has_count_files = True
            has_fastq = False  # Skip FASTQ processing if we have a count file
        except Exception as e:
            logging.error(f"Error copying count file: {e}")
            return {"error": f"Failed to copy count file: {str(e)}"}
    else:
        # Determine input data types - look directly in experiment subdirectory or data dir
        progress_reporter.update("File Discovery", "Looking for FASTQ and count files")
        # Check for FASTQ files directly in the input directory
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
        
        # Determine if we're processing FASTQ files
        has_fastq = len(fastq_files) > 0
        
        # Check for FASTQ files to determine if library file is required
        has_fastq = False
        
        # Only check for FASTQ files if no count file is provided
        if not count_file:
            # Check for FASTQ files directly in the input directory
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
            
            has_fastq = len(fastq_files) > 0
        
        # Check if library file exists - only required for FASTQ processing
        if has_fastq and not os.path.exists(library_file):
            # Also check legacy location
            legacy_library = os.path.join(input_dir, "library.csv")
            if os.path.exists(legacy_library):
                library_file = legacy_library
                progress.update("Library File", f"Using legacy library file: {library_file}")
            else:
                print(f"Error: Library file not found at {library_file} or {legacy_library}")
                print("A library file is required when processing FASTQ files.")
                return 1
        elif has_fastq:
            progress.update("Library File", f"Using library file: {library_file}")
        else:
            # For count files, library file is optional
            if os.path.exists(library_file):
                progress.update("Library File", f"Using optional library file: {library_file}")
            else:
                progress.update("Library File", "No library file found (not required for count files)")
        
        if args.contrasts_file is None and not os.path.exists(contrasts_file):
            # Also check legacy location
            legacy_contrasts = os.path.join(input_dir, "contrasts.csv")
            if os.path.exists(legacy_contrasts):
                contrasts_file = legacy_contrasts
                progress.update("Contrast File", f"Using legacy contrasts file: {contrasts_file}")
            else:
                print(f"Warning: No contrasts file found at {contrasts_file} or {legacy_contrasts}")
                print("Differential analysis will be skipped.")
        
        # Set logging level
        log_level = logging.DEBUG if args.verbose else logging.INFO
        logging.getLogger().setLevel(log_level)
        
        # Determine parallel processing
        # Check if running under Snakemake and disable parallelism if so
        running_in_snakemake = 'SNAKEMAKE' in os.environ or 'snakemake' in os.environ
        if running_in_snakemake:
            logging.info("Detected Snakemake execution environment - disabling internal parallelism")
            parallel = False
        else:
            parallel = not args.no_parallel
        
        # Run the pipeline
        progress.update("Pipeline Start", "Running main pipeline")
        results = run_pipeline(
            input_dir=input_dir,
            output_dir=output_dir,
            library_file=library_file,
            experiment_name=args.experiment_name,
            contrasts_file=contrasts_file,
            sample_sheet=args.sample_sheet,
            norm_method=args.norm_method,
            fdr_threshold=args.fdr_threshold,
            overwrite=args.overwrite,
            skip_drugz=args.skip_drugz,
            skip_qc=args.skip_qc,
            skip_mle=args.skip_mle,
            use_docker=True,  # Always use Docker
            count_file=count_file,
            parallel=parallel,
            max_workers=args.workers,
            progress_reporter=progress
        )
        
        if "error" in results:
            print(f"Error: {results['error']}")
            return 1
        
        # Final progress update
        progress.update("Completion", "Pipeline completed successfully")
        
        print(f"Results directory: {output_dir}")
        
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
    parser.add_argument("--count-file", help="Path to a count file (.count or .csv) to use directly instead of looking in input-dir")
    
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
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    parser.add_argument("--no-parallel", action="store_true", help="Disable parallel processing")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (default: 4)")
    
    args = parser.parse_args()
    
    # Create absolute paths
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    # Import progress reporter
    from analysis_pipeline.core.logging_setup import ProgressReporter
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a progress reporter
    # Estimate the total number of steps based on arguments
    total_steps = 5  # Base steps (initialization, finding files, validation, etc.)
    if args.count_file:
        total_steps += 1  # Processing count file
    else:
        total_steps += 3  # Finding FASTQ/count files + processing + merging
    
    # Add steps for each contrast (estimated if contrasts file is provided)
    if args.contrasts_file and os.path.exists(args.contrasts_file):
        try:
            contrasts_df = pd.read_csv(args.contrasts_file)
            if "contrast" in contrasts_df.columns:
                num_contrasts = len(contrasts_df["contrast"])
                # Each contrast has RRA, possibly MLE, possibly DrugZ
                contrast_steps = num_contrasts
                if not args.skip_mle:
                    contrast_steps += num_contrasts
                if not args.skip_drugz:
                    contrast_steps += num_contrasts
                total_steps += contrast_steps
        except:
            total_steps += 5  # Fallback estimate for contrasts
    
    # Create progress reporter
    progress = ProgressReporter(total_steps, args.experiment_name)
    progress.update("Initialization", "Setting up pipeline")
    
    # Default library and contrasts files within the experiment directory
    library_file = args.library_file or os.path.join(input_dir, args.experiment_name, "library.csv")
    contrasts_file = args.contrasts_file or os.path.join(input_dir, args.experiment_name, "contrasts.csv")
    
    progress.update("Input Validation", "Checking input files")
    
    # Handle count file if provided
    count_file = None
    if args.count_file:
        if not os.path.exists(args.count_file):
            print(f"Error: Count file not found at {args.count_file}")
            return 1
        
        count_file = os.path.abspath(args.count_file)
        progress.update("Count File", f"Using provided count file: {count_file}")
        
        # If it's a CSV file, convert it to tab-delimited format
        if count_file.endswith('.csv'):
            try:
                from analysis_pipeline.core.file_handling import ensure_count_formats
                progress.update("CSV Conversion", f"Converting CSV to tab-delimited: {count_file}")
                csv_path, tab_path = ensure_count_formats(count_file, output_dir)
                count_file = tab_path  # Use the tab-delimited file for analysis
                print(f"Created both formats - CSV: {csv_path}, Tab-delimited: {tab_path}")
            except Exception as e:
                print(f"Error converting count file formats: {str(e)}")
                return 1
        # If it's a tab-delimited file, ensure we also have a CSV version
        elif count_file.endswith('.count') or count_file.endswith('.counts') or count_file.endswith('.txt'):
            try:
                from analysis_pipeline.core.file_handling import ensure_count_formats
                progress.update("Format Conversion", f"Ensuring both formats exist for: {count_file}")
                csv_path, tab_path = ensure_count_formats(count_file, output_dir)
                count_file = tab_path  # Continue using the tab-delimited file for analysis
                print(f"Created both formats - CSV: {csv_path}, Tab-delimited: {tab_path}")
            except Exception as e:
                print(f"Error ensuring count file formats: {str(e)}")
                # Continue with the original file
    
    # Check for FASTQ files to determine if library file is required
    has_fastq = False
    
    # Only check for FASTQ files if no count file is provided
    if not count_file:
        # Check for FASTQ files directly in the input directory
        fastq_pattern = os.path.join(input_dir, args.experiment_name, "**", "*.fastq*")
        fastq_files = glob.glob(fastq_pattern, recursive=True)
        
        # If no fastq files found in experiment dir, check in data dir
        if not fastq_files:
            fastq_pattern = os.path.join(input_dir, args.experiment_name, "data", "**", "*.fastq*")
            fastq_files = glob.glob(fastq_pattern, recursive=True)
        
        # Also look for fastq in fastq subdirectory
        if not fastq_files:
            fastq_pattern = os.path.join(input_dir, args.experiment_name, "fastq", "**", "*.fastq*")
            fastq_files = glob.glob(fastq_pattern, recursive=True)
        
        # Also try fq extension
        if not fastq_files:
            fastq_pattern = os.path.join(input_dir, args.experiment_name, "**", "*.fq*")
            fastq_files = glob.glob(fastq_pattern, recursive=True)
        
        has_fastq = len(fastq_files) > 0
    
    # Check if library file exists - only required for FASTQ processing
    if has_fastq and not os.path.exists(library_file):
        # Also check legacy location
        legacy_library = os.path.join(input_dir, "library.csv")
        if os.path.exists(legacy_library):
            library_file = legacy_library
            progress.update("Library File", f"Using legacy library file: {library_file}")
        else:
            print(f"Error: Library file not found at {library_file} or {legacy_library}")
            print("A library file is required when processing FASTQ files.")
            return 1
    elif has_fastq:
        progress.update("Library File", f"Using library file: {library_file}")
    else:
        # For count files, library file is optional
        if os.path.exists(library_file):
            progress.update("Library File", f"Using optional library file: {library_file}")
        else:
            progress.update("Library File", "No library file found (not required for count files)")
    
    if args.contrasts_file is None and not os.path.exists(contrasts_file):
        # Also check legacy location
        legacy_contrasts = os.path.join(input_dir, "contrasts.csv")
        if os.path.exists(legacy_contrasts):
            contrasts_file = legacy_contrasts
            progress.update("Contrast File", f"Using legacy contrasts file: {contrasts_file}")
        else:
            print(f"Warning: No contrasts file found at {contrasts_file} or {legacy_contrasts}")
            print("Differential analysis will be skipped.")
    
    # Set logging level
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(log_level)
    
    # Determine parallel processing
    # Check if running under Snakemake and disable parallelism if so
    running_in_snakemake = 'SNAKEMAKE' in os.environ or 'snakemake' in os.environ
    if running_in_snakemake:
        logging.info("Detected Snakemake execution environment - disabling internal parallelism")
        parallel = False
    else:
        parallel = not args.no_parallel
    
    # Run the pipeline
    progress.update("Pipeline Start", "Running main pipeline")
    results = run_pipeline(
        input_dir=input_dir,
        output_dir=output_dir,
        library_file=library_file,
        experiment_name=args.experiment_name,
        contrasts_file=contrasts_file,
        sample_sheet=args.sample_sheet,
        norm_method=args.norm_method,
        fdr_threshold=args.fdr_threshold,
        overwrite=args.overwrite,
        skip_drugz=args.skip_drugz,
        skip_qc=args.skip_qc,
        skip_mle=args.skip_mle,
        use_docker=True,  # Always use Docker
        count_file=count_file,
        parallel=parallel,
        max_workers=args.workers,
        progress_reporter=progress
    )
    
    if "error" in results:
        print(f"Error: {results['error']}")
        return 1
    
    # Final progress update
    progress.update("Completion", "Pipeline completed successfully")
    
    print(f"Results directory: {output_dir}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 