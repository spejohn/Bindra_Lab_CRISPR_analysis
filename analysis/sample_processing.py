"""
Sample processing functions for CRISPR analysis pipeline.
Handles individual sample processing including count generation from FASTQ files.
"""

import os
import logging
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any

# Import from core modules
from analysis_pipeline.core.file_handling import ensure_output_dir, reverse_dir
from analysis_pipeline.docker.docker_utils import run_docker_container, verify_docker
from analysis_pipeline.core.config import DOCKER_IMAGE, FASTQ_PATTERNS, SAMPLE_SHEET_REQUIRED_COLUMNS


def mageck_count_single_sample(
    fastq_file: str,
    library_file: str,
    output_dir: str,
    sample_name: Optional[str] = None,
    count_options: Optional[Dict[str, Any]] = None
) -> Tuple[bool, str]:
    """
    Run MAGeCK count on a single FASTQ file using Docker.
    
    Args:
        fastq_file: Path to the FASTQ file
        library_file: Path to the library file
        output_dir: Directory for output files
        sample_name: Sample name (defaults to fastq filename without extension)
        count_options: Additional options for MAGeCK count
        
    Returns:
        Tuple of (success, output_file_path)
    """
    # Verify Docker is available
    if not verify_docker():
        error_msg = "Docker verification failed, cannot run MAGeCK count"
        logging.error(error_msg)
        return (False, error_msg)
        
    # Ensure the output directory exists
    ensure_output_dir(output_dir)
    
    # Set default count options
    if count_options is None:
        count_options = {}
    
    # Determine sample name if not provided
    if sample_name is None:
        sample_name = Path(fastq_file).stem
        # Remove .fastq or .fq if present in the stem
        for suffix in [".fastq", ".fq"]:
            if sample_name.endswith(suffix):
                sample_name = sample_name[:-len(suffix)]
        # Remove .gz if present
        if sample_name.endswith(".gz"):
            sample_name = sample_name[:-3]
            
    logging.info(f"Processing sample: {sample_name} from {fastq_file}")
    
    # Create output filename
    output_file = os.path.join(output_dir, f"{sample_name}.count.txt")
    
    # Prepare Docker volumes
    fastq_dir = os.path.dirname(os.path.abspath(fastq_file))
    library_dir = os.path.dirname(os.path.abspath(library_file))
    
    volumes = {
        fastq_dir: {"bind": "/fastq", "mode": "ro"},
        library_dir: {"bind": "/library", "mode": "ro"},
        output_dir: {"bind": "/output", "mode": "rw"}
    }
    
    # Convert paths to Docker paths
    docker_fastq = f"/fastq/{os.path.basename(fastq_file)}"
    docker_library = f"/library/{os.path.basename(library_file)}"
    docker_output = f"/output/{os.path.basename(output_file)}"
    
    # Prepare MAGeCK count command
    command = ["mageck", "count"]
    command.extend(["--fastq", docker_fastq])
    command.extend(["--list-seq", docker_library])
    command.extend(["--sample-label", sample_name])
    command.extend(["--output-prefix", f"/output/{sample_name}"])
    
    # Add additional options
    for option, value in count_options.items():
        if value is None:
            command.append(f"--{option}")
        else:
            command.extend([f"--{option}", str(value)])
    
    # Run Docker container
    exit_code, output = run_docker_container(
        command=command,
        volumes=volumes,
        stream_logs=True
    )
    
    if exit_code != 0:
        logging.error(f"MAGeCK count failed for sample {sample_name}")
        return (False, f"MAGeCK count failed: {output}")
    
    # Check if output file exists
    if not os.path.exists(output_file):
        error_msg = f"Expected output file {output_file} not found after MAGeCK count"
        logging.error(error_msg)
        return (False, error_msg)
    
    count_log_file = os.path.join(output_dir, f"{sample_name}.count.log")
    if os.path.exists(count_log_file):
        logging.info(f"MAGeCK count log file created: {count_log_file}")
    
    logging.info(f"MAGeCK count completed successfully for sample {sample_name}")
    logging.info(f"Output count file: {output_file}")
    
    return (True, output_file)


def generate_sample_sheet(
    fastq_dir: str,
    output_dir: str,
    contrast_name: str = "contrast",
    include_patterns: Optional[List[str]] = None,
    experiment_name: Optional[str] = None
) -> Optional[pd.DataFrame]:
    """
    Generate a sample sheet based on FASTQ files in the input directory.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        output_dir: Directory to save the sample sheet
        contrast_name: Name prefix for the sample sheet file
        include_patterns: File patterns to include (defaults to FASTQ_PATTERNS from config)
        experiment_name: Name of the experiment to look for in a subdirectory
        
    Returns:
        DataFrame containing the sample sheet or None if no files found
    """
    if include_patterns is None:
        include_patterns = FASTQ_PATTERNS
    
    # If experiment_name is provided, check in experiment subdirectory
    if experiment_name:
        experiment_path = Path(fastq_dir) / experiment_name
        if experiment_path.exists():
            logging.info(f"Looking for FASTQ files in experiment directory: {experiment_path}")
            
            # Check for a data directory within the experiment directory
            data_dir = experiment_path / "data"
            if data_dir.exists() and data_dir.is_dir():
                logging.info(f"Found data directory: {data_dir}")
                fastq_dir = str(data_dir)
            else:
                # Also check for a fastq directory
                fastq_subdir = experiment_path / "fastq"
                if fastq_subdir.exists() and fastq_subdir.is_dir():
                    logging.info(f"Found fastq directory: {fastq_subdir}")
                    fastq_dir = str(fastq_subdir)
                else:
                    # Use the experiment directory directly
                    fastq_dir = str(experiment_path)
    
    logging.info(f"Generating sample sheet from FASTQ files in {fastq_dir}")
    logging.info(f"Using file patterns: {include_patterns}")
    
    # Find all FASTQ files
    fastq_files = []
    for root, _, files in os.walk(fastq_dir):
        for file in files:
            for pattern in include_patterns:
                if file.lower().endswith(pattern.lower()):
                    fastq_files.append(os.path.join(root, file))
                    break
    
    if not fastq_files:
        logging.warning(f"No FASTQ files found in {fastq_dir}")
        return None
    
    logging.info(f"Found {len(fastq_files)} FASTQ files")
    
    # Create sample sheet
    samples = []
    for fastq_file in fastq_files:
        sample_name = Path(fastq_file).stem
        # Remove .fastq or .fq if present in the stem
        for suffix in [".fastq", ".fq"]:
            if sample_name.lower().endswith(suffix):
                sample_name = sample_name[:-len(suffix)]
        
        # Remove .gz if present
        if sample_name.lower().endswith(".gz"):
            sample_name = sample_name[:-3]
            
        rel_path = os.path.relpath(fastq_file, fastq_dir)
        
        samples.append({
            "sample": sample_name,
            "fastq_path": rel_path,
            "condition": "unknown",  # Default condition, to be filled by user
            "replicate": 1           # Default replicate, to be filled by user
        })
    
    # Create DataFrame
    sample_sheet = pd.DataFrame(samples)
    
    # Ensure all required columns are present
    for col in SAMPLE_SHEET_REQUIRED_COLUMNS:
        if col not in sample_sheet.columns:
            sample_sheet[col] = ""
    
    # Save to file using the contrast name in the filename
    sample_sheet_path = os.path.join(output_dir, f"{contrast_name}_samples.txt")
    sample_sheet.to_csv(sample_sheet_path, index=False)
    
    logging.info(f"Generated sample sheet with {len(samples)} samples")
    logging.info(f"Sample sheet saved to {sample_sheet_path}")
    
    return sample_sheet


def process_all_samples(
    fastq_dir: str,
    library_file: str,
    output_dir: str,
    sample_sheet: Optional[pd.DataFrame] = None,
    count_options: Optional[Dict[str, Any]] = None,
    overwrite: bool = False,
    experiment_name: Optional[str] = None,
    parallel: bool = True,
    max_workers: int = 4
) -> Dict[str, str]:
    """
    Process all FASTQ samples in a directory to count files.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        library_file: Path to the library file
        output_dir: Directory for output files
        sample_sheet: DataFrame containing sample information (will be generated if None)
        count_options: Options for the count function
        overwrite: Whether to overwrite existing count files
        experiment_name: Name of the experiment to look for in a subdirectory
        parallel: Whether to use parallel processing
        max_workers: Maximum number of parallel workers
        
    Returns:
        Dictionary mapping sample names to count file paths
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # If sample_sheet is not provided, generate one
    if sample_sheet is None:
        # Generate sample sheet for all FASTQ files
        experiment_base = os.path.basename(output_dir)
        sample_sheet = generate_sample_sheet(fastq_dir, output_dir, experiment_base, experiment_name=experiment_name)
        
        if sample_sheet is None:
            logging.error(f"No FASTQ files found in {fastq_dir}")
            return {}
    
    # Check if fastq_dir exists
    if not os.path.exists(fastq_dir):
        logging.error(f"FASTQ directory {fastq_dir} does not exist")
        return {}
    
    # Adjust fastq_dir if experiment_name is provided
    if experiment_name:
        # First, check if the experiment directory exists
        experiment_path = Path(fastq_dir) / experiment_name
        if experiment_path.exists():
            # Check for data directory
            data_dir = experiment_path / "data"
            if data_dir.exists() and data_dir.is_dir():
                fastq_dir = str(data_dir)
            else:
                # Check for fastq directory
                fastq_subdir = experiment_path / "fastq"
                if fastq_subdir.exists() and fastq_subdir.is_dir():
                    fastq_dir = str(fastq_subdir)
                else:
                    # Use experiment directory as is
                    fastq_dir = str(experiment_path)
                    
            logging.info(f"Using experiment-specific directory for FASTQ files: {fastq_dir}")
        else:
            logging.warning(f"Experiment directory {experiment_path} not found")
    
    # Count files dict to return
    count_files = {}
    
    # Check if we should process in parallel
    if parallel and len(sample_sheet) > 1:
        import concurrent.futures
        # Determine how many workers to use
        n_workers = min(len(sample_sheet), max_workers)
        logging.info(f"Processing {len(sample_sheet)} samples in parallel with {n_workers} workers")
        
        # Function to process a single sample in parallel
        def process_single_sample(row_tuple):
            i, row = row_tuple
            sample_name = row["sample"]
            
            if "fastq_path" not in row:
                logging.warning(f"Missing fastq_path for sample {sample_name}, skipping")
                return (sample_name, None)
                
            fastq_path = os.path.join(fastq_dir, row["fastq_path"])
            count_file = os.path.join(output_dir, f"{sample_name}.count")
            
            logging.info(f"Processing sample {sample_name} ({i+1}/{len(sample_sheet)})")
            
            # Check if count file already exists and whether to overwrite
            if os.path.exists(count_file) and not overwrite:
                logging.info(f"Count file {count_file} already exists, skipping (use --overwrite to force recount)")
                return (sample_name, count_file)
            
            # Check if FASTQ file exists
            if not os.path.exists(fastq_path):
                logging.error(f"FASTQ file {fastq_path} does not exist")
                return (sample_name, None)
            
            # Count FASTQ file
            success, result_file = mageck_count_single_sample(
                fastq_file=fastq_path,
                library_file=library_file,
                output_dir=output_dir,
                sample_name=sample_name,
                count_options=count_options
            )
            
            if success:
                logging.info(f"Counted {sample_name} -> {result_file}")
                return (sample_name, result_file)
            else:
                logging.error(f"Failed to count {sample_name}")
                return (sample_name, None)
        
        # Process samples in parallel using ThreadPoolExecutor
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {}
            
            # Submit all tasks
            for i, row in sample_sheet.iterrows():
                future = executor.submit(process_single_sample, (i, row))
                futures[future] = row["sample"]
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(futures):
                sample_name, result_file = future.result()
                if result_file:
                    count_files[sample_name] = result_file
    else:
        # Process each sample sequentially
        for i, row in sample_sheet.iterrows():
            sample_name = row["sample"]
            if "fastq_path" not in row:
                logging.warning(f"Missing fastq_path for sample {sample_name}, skipping")
                continue
                
            fastq_path = os.path.join(fastq_dir, row["fastq_path"])
            count_file = os.path.join(output_dir, f"{sample_name}.count")
            
            logging.info(f"Processing sample {sample_name} ({i+1}/{len(sample_sheet)})")
            
            # Check if count file already exists and whether to overwrite
            if os.path.exists(count_file) and not overwrite:
                logging.info(f"Count file {count_file} already exists, skipping (use --overwrite to force recount)")
                count_files[sample_name] = count_file
                continue
            
            # Check if FASTQ file exists
            if not os.path.exists(fastq_path):
                logging.error(f"FASTQ file {fastq_path} does not exist")
                continue
            
            # Count FASTQ file
            success, result_file = mageck_count_single_sample(
                fastq_file=fastq_path,
                library_file=library_file,
                output_dir=output_dir,
                sample_name=sample_name,
                count_options=count_options
            )
            
            if success:
                count_files[sample_name] = result_file
                logging.info(f"Counted {sample_name} -> {result_file}")
            else:
                logging.error(f"Failed to count {sample_name}")
    
    logging.info(f"Processed {len(count_files)}/{len(sample_sheet)} samples")
    
    return count_files


def merge_count_files(
    count_files: Dict[str, str],
    output_dir: str,
    contrast_name: str = "contrast"
) -> Optional[str]:
    """
    Merge multiple count files into a single count table.
    
    Args:
        count_files: Dictionary mapping sample names to count file paths
        output_dir: Directory to save the merged file
        contrast_name: Name of the contrast (used for file naming)
        
    Returns:
        Path to the merged count file or None if failed
    """
    if not count_files:
        logging.error("No count files provided for merging")
        return None
    
    logging.info(f"Merging {len(count_files)} count files")
    
    try:
        # Process each count file to extract counts
        merged_data = {}
        
        for sample, file_path in count_files.items():
            logging.info(f"Processing count file: {file_path}")
            
            # Read the count file
            try:
                df = pd.read_csv(file_path, sep='\t')
                
                # Extract sgRNA IDs and counts
                merged_data[sample] = dict(zip(df['sgRNA'], df['count']))
                
                logging.info(f"Added {len(df)} guides from sample {sample}")
                
            except Exception as e:
                logging.error(f"Error reading count file {file_path}: {e}")
                continue
        
        if not merged_data:
            logging.error("No valid count data found in any files")
            return None
        
        # Get all unique sgRNAs
        all_sgrnas = set()
        for sample_data in merged_data.values():
            all_sgrnas.update(sample_data.keys())
        
        logging.info(f"Found {len(all_sgrnas)} unique sgRNAs across all samples")
        
        # Create a DataFrame with all sgRNAs and counts
        data = {'sgRNA': sorted(all_sgrnas)}
        
        for sample in sorted(merged_data.keys()):
            sample_counts = merged_data[sample]
            data[sample] = [sample_counts.get(sgrna, 0) for sgrna in data['sgRNA']]
        
        # Create merged DataFrame
        merged_df = pd.DataFrame(data)
        
        # Save the merged count file directly in the experiment directory
        output_file = os.path.join(output_dir, f"{contrast_name}.count")
        merged_df.to_csv(output_file, sep='\t', index=False)
        
        logging.info(f"Merged count file saved to {output_file}")
        
        return output_file
        
    except Exception as e:
        logging.error(f"Error merging count files: {e}")
        return None 