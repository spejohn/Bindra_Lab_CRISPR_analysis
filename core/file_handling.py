"""
File and directory handling utilities for CRISPR analysis pipeline.
"""

import os
import re
import logging
import shutil
import pandas as pd
from pathlib import Path
from typing import Set, List, Dict, Any, Optional, Union

from analysis_pipeline.core.config import (
    CONTRAST_TABLE_PATTERN,
    READ_COUNT_PATTERN,
    FASTQ_PATTERNS
)


def reverse_dir(input_dir: str, root: str, output_dir: str) -> Path:
    """
    Create an output directory structure mirroring the input directory structure.
    
    Args:
        input_dir: Base input directory
        root: Current directory being processed
        output_dir: Base output directory
        
    Returns:
        Path to the mirrored directory in the output directory
    """
    try:
        # Use regex to find and remove the input_dir path from root
        if input_dir in root:
            dir_str = re.sub(f'^{re.escape(input_dir)}', '', root)
        else:
            logging.warning(f"Input directory string not found for {root} file.")
            dir_str = Path(root).name

        # Combine the paths to mirror the input directory
        out_path = Path(output_dir) / Path(dir_str.lstrip('\\/'))

        # Create dir if doesn't exist
        os.makedirs(out_path, exist_ok=True)
        
        return out_path
        
    except Exception as e:
        logging.error(f"Error creating mirrored directory: {str(e)}")
        # Fallback to using just the name of the root directory
        fallback_path = Path(output_dir) / Path(root).name
        os.makedirs(fallback_path, exist_ok=True)
        return fallback_path


def find_input_files(input_dir: str, experiment_name: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
    """
    Scan the input directory for FASTQ files and read count tables.
    
    Args:
        input_dir: Input directory to scan
        experiment_name: Optional experiment name to look for in a specific subdirectory
        
    Returns:
        Dictionary with input files organized by type
    """
    input_path = Path(input_dir)
    
    # If experiment_name is provided, focus on that subdirectory
    if experiment_name:
        experiment_path = input_path / experiment_name
        if experiment_path.exists():
            input_path = experiment_path
            logging.info(f"Using experiment-specific directory: {input_path}")
        else:
            logging.warning(f"Experiment directory {experiment_path} not found. Using {input_path} instead.")
    
    # Dictionary to store input files
    input_files = {
        'fastq_dirs': {},  # Format: {dir_path: {'contrast_table': path, 'design_matrix': path, 'fastq_files': [paths]}}
        'read_counts': {},  # Format: {file_path: {'contrast_table': path, 'design_matrix': path}}
        'library_file': None,
        'contrast_file': None,
        'design_matrix': None
    }
    
    # Define patterns for recognizing design matrix files
    DESIGN_MATRIX_PATTERN = "*design*.txt"
    LIBRARY_PATTERNS = ["*library*.csv", "*library*.txt", "*guide*.csv", "*guide*.txt"]
    
    fastq_dirs = {}  # Track FASTQ directories we've already processed
    
    # First look for library, contrast, and design matrix files at the experiment level
    for pattern in LIBRARY_PATTERNS:
        library_files = list(input_path.glob(pattern))
        if library_files:
            input_files['library_file'] = str(library_files[0])
            logging.info(f"Found library file: {input_files['library_file']}")
            break
    
    contrast_files = list(input_path.glob(CONTRAST_TABLE_PATTERN))
    if contrast_files:
        input_files['contrast_file'] = str(contrast_files[0])
        logging.info(f"Found contrast file: {input_files['contrast_file']}")
    
    design_files = list(input_path.glob(DESIGN_MATRIX_PATTERN))
    if design_files:
        input_files['design_matrix'] = str(design_files[0])
        logging.info(f"Found design matrix: {input_files['design_matrix']}")
        
    # We only recognize "fastq" and "counts" directories
    # Do NOT look for generic "data" directory as it causes ambiguity

    # Walk through the input directory
    for root, dirs, files in os.walk(input_path):
        root_path = Path(root)
        
        # Skip if no files in directory
        if not files:
            continue
        
        # Find the contrast table (if any) in this directory
        cnttbl_files = [f for f in files if fnmatch_lower(f, CONTRAST_TABLE_PATTERN)]
        cnttbl_path = str(root_path / cnttbl_files[0]) if cnttbl_files else input_files.get('contrast_file')
        
        # Find the design matrix (if any) in this directory
        design_files = [f for f in files if fnmatch_lower(f, DESIGN_MATRIX_PATTERN)]
        design_path = str(root_path / design_files[0]) if design_files else input_files.get('design_matrix')
        
        # Check if this is a read count directory
        if "read_count" in root.lower() or any(f.lower().endswith('.count') for f in files):
            # Find read count files (RC files)
            rc_files = [f for f in files if fnmatch_lower(f, READ_COUNT_PATTERN)]
            
            if rc_files and cnttbl_path:
                # Process each read count file
                for rc_file in rc_files:
                    rc_path = str(root_path / rc_file)
                    input_files['read_counts'][rc_path] = {
                        'contrast_table': cnttbl_path,
                        'design_matrix': design_path
                    }
                
                logging.info(f"Found {len(rc_files)} read count files in {root}")
                
        # Check if this is a FASTQ directory
        elif "fastq" in root.lower() and root not in fastq_dirs:
            # Check for FASTQ files
            fastq_files = []
            for pattern in FASTQ_PATTERNS:
                fastq_files.extend([str(root_path / f) for f in files if fnmatch_lower(f, pattern)])
            
            if fastq_files:
                fastq_dirs[root] = True  # Mark this directory as processed
                input_files['fastq_dirs'][root] = {
                    'contrast_table': cnttbl_path,
                    'design_matrix': design_path,
                    'fastq_files': fastq_files
                }
                
                logging.info(f"Found {len(fastq_files)} FASTQ files in {root}")
    
    # Count the number of design matrices found
    design_matrices_count = sum(1 for _, info in input_files['read_counts'].items() if info.get('design_matrix') is not None)
    design_matrices_count += sum(1 for _, info in input_files['fastq_dirs'].items() if info.get('design_matrix') is not None)
    
    logging.info(f"Found {len(input_files['read_counts'])} read count files and {len(input_files['fastq_dirs'])} FASTQ directories")
    logging.info(f"Found {design_matrices_count} design matrices for MLE analysis")
    
    return input_files


def fnmatch_lower(filename: str, pattern: str) -> bool:
    """
    Case-insensitive file pattern matching.
    
    Args:
        filename: Filename to match
        pattern: Pattern to match against (with wildcards)
        
    Returns:
        True if the filename matches the pattern, False otherwise
    """
    import fnmatch
    return fnmatch.fnmatch(filename.lower(), pattern.lower())


def ensure_output_dir(output_dir: str) -> None:
    """
    Ensure the output directory exists.
    
    Args:
        output_dir: Path to the output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Output directory created/verified: {output_dir}")


def make_count_table(rc_file: str, output_dir: Optional[str] = None) -> str:
    """
    Convert a read count table from CSV to tab-delimited format for MAGeCK.
    
    Args:
        rc_file: Path to read count CSV file
        output_dir: Optional output directory (defaults to same as rc_file)
        
    Returns:
        Path to the tab-delimited count table
    """
    try:
        rc_path = Path(rc_file)
        out_dir = Path(output_dir) if output_dir else rc_path.parent
        
        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)
        
        # Load the read count table
        rc_data = pd.read_csv(rc_path)
        logging.info(f"Loaded read count table: {rc_file} with {len(rc_data)} rows")
        
        # Create the count table filename
        count_file = out_dir / rc_path.name.replace(".csv", ".count.txt")
        
        # Save as tab-delimited
        rc_data.to_csv(count_file, sep="\t", index=False)
        logging.info(f"Created tab-delimited count table: {count_file}")
        
        return str(count_file)
        
    except Exception as e:
        logging.error(f"Error creating count table from {rc_file}: {str(e)}")
        raise


def create_experiment_dirs(base_dir: str, contrast_name: str) -> Dict[str, str]:
    """
    Create standard directory structure for a contrast.
    
    Args:
        base_dir: Base directory for the contrast
        contrast_name: Name of the contrast
        
    Returns:
        Dictionary with paths to created directories
    """
    base_path = Path(base_dir)
    
    # Simplified structure - only use the contrast directory itself
    # without creating subdirectories for each analysis type
    dirs = {
        'contrast': str(base_path)  # The contrast directory serves as the root for all files
    }
    
    # Create the directory
    os.makedirs(base_path, exist_ok=True)
    
    logging.info(f"Created directory structure for contrast: {contrast_name}")
    return dirs


def identify_analyzed_experiments(output_dir: str) -> Dict[str, Dict[str, bool]]:
    """
    Identify experiments that have already been analyzed in the output directory.
    
    Args:
        output_dir: Path to the output directory
        
    Returns:
        Dictionary mapping experiment names to status of analyses
    """
    analyzed_experiments = {}
    
    if not os.path.exists(output_dir):
        return analyzed_experiments
    
    # Check for experiment directories
    for exp_dir in os.listdir(output_dir):
        exp_path = os.path.join(output_dir, exp_dir)
        
        if not os.path.isdir(exp_path):
            continue
        
        # Initialize status for this experiment
        analyzed_experiments[exp_dir] = {
            "count_files": False,
            "rra_analysis": False,
            "mle_analysis": False,
            "drugz_analysis": False,
            "qc_plots": False
        }
        
        # Check for count files
        count_dir = os.path.join(exp_path, "counts")
        if os.path.exists(count_dir) and os.listdir(count_dir):
            analyzed_experiments[exp_dir]["count_files"] = True
        
        # Check for RRA analysis
        rra_dir = os.path.join(exp_path, "RRA")
        if os.path.exists(rra_dir) and os.listdir(rra_dir):
            analyzed_experiments[exp_dir]["rra_analysis"] = True
        
        # Check for MLE analysis
        mle_dir = os.path.join(exp_path, "MLE")
        if os.path.exists(mle_dir) and os.listdir(mle_dir):
            analyzed_experiments[exp_dir]["mle_analysis"] = True
        
        # Check for DrugZ analysis
        drugz_dir = os.path.join(exp_path, "drugz")
        if os.path.exists(drugz_dir) and os.listdir(drugz_dir):
            analyzed_experiments[exp_dir]["drugz_analysis"] = True
        
        # Check for QC plots
        qc_dir = os.path.join(exp_path, "qc")
        if os.path.exists(qc_dir) and os.listdir(qc_dir):
            analyzed_experiments[exp_dir]["qc_plots"] = True
    
    return analyzed_experiments


def convert_results_to_csv(result_file: str, analysis_type: str) -> str:
    """
    Convert analysis results (RRA or DrugZ) to CSV format with appropriate suffix.
    
    Args:
        result_file: Path to the result file (tab-delimited)
        analysis_type: Type of analysis ('RRA' or 'DrugZ')
        
    Returns:
        Path to the created CSV file
    """
    try:
        # Determine suffix based on analysis type
        if analysis_type.upper() == 'RRA':
            suffix = '_gMGK'
        elif analysis_type.upper() == 'DRUGZ':
            suffix = '_gDZ'
        else:
            suffix = f'_{analysis_type}'
        
        # Get the base name and directory
        result_path = Path(result_file)
        result_dir = result_path.parent
        
        # Create new filename with appropriate suffix
        # Remove old suffix if present
        base_name = result_path.stem
        for old_suffix in ['_RRA', '_DrugZ', '.gene_summary']:
            if old_suffix in base_name:
                base_name = base_name.replace(old_suffix, '')
        
        csv_file = os.path.join(result_dir, f"{base_name}{suffix}.csv")
        
        # Read the tab-delimited file and convert to CSV
        df = pd.read_csv(result_file, sep='\t')
        df.to_csv(csv_file, index=False)
        
        logging.info(f"Converted {result_file} to CSV format: {csv_file}")
        
        return csv_file
    
    except Exception as e:
        logging.error(f"Error converting {result_file} to CSV: {str(e)}")
        return result_file 