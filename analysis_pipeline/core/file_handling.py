"""
File and directory handling utilities for CRISPR analysis pipeline.
"""

import os
import re
import logging
import shutil
import pandas as pd
from pathlib import Path
from typing import Set, List, Dict, Any, Optional, Union, Tuple

from analysis_pipeline.core.config import (
    CONTRAST_TABLE_PATTERN,
    READ_COUNT_PATTERN,
    FASTQ_PATTERNS,
    COUNT_CSV_PATTERN,
    CONTRAST_REQUIRED_COLUMNS
)
from analysis_pipeline.core.utils import retry_operation

# Define patterns
DESIGN_MATRIX_PATTERN = "*.txt"


@retry_operation(max_attempts=3, delay=2)
def make_count_table(rc_file: str, output_dir: Optional[str] = None) -> str:
    """
    Convert a read count table from CSV to tab-delimited format for MAGeCK.
    Uses retry_operation decorator to handle temporary file access issues.
    
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


def fnmatch_lower(name: str, pattern: str) -> bool:
    """
    Case-insensitive version of fnmatch.
    
    Args:
        name: String to match
        pattern: Pattern to match against
        
    Returns:
        True if name matches pattern, False otherwise
    """
    import fnmatch
    return fnmatch.fnmatch(name.lower(), pattern.lower())


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


@retry_operation(max_attempts=3, delay=2)
def ensure_count_formats(count_file_path: str, output_dir: Optional[str] = None) -> Tuple[str, str]:
    """
    Ensures that a count file exists in both CSV and tab-delimited formats.
    If given a CSV, creates a tab-delimited version.
    If given a tab-delimited file, creates a CSV version.
    
    Args:
        count_file_path: Path to the count file (either CSV or tab-delimited)
        output_dir: Optional output directory (defaults to same as count_file_path)
        
    Returns:
        Tuple of (csv_file_path, tab_file_path)
    """
    try:
        file_path = Path(count_file_path)
        out_dir = Path(output_dir) if output_dir else file_path.parent
        
        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)
        
        # Determine the file type and create the other format
        if file_path.name.lower().endswith('.csv'):
            # This is a CSV file - create tab-delimited version
            csv_file_path = str(file_path)
            # Load the CSV file
            data = pd.read_csv(file_path)
            
            # Create the tab-delimited filename
            if '.counts.csv' in file_path.name.lower():
                tab_file_path = str(out_dir / file_path.name.replace('.counts.csv', '.counts.txt'))
            elif '.fastq_counts.csv' in file_path.name.lower():
                tab_file_path = str(out_dir / file_path.name.replace('.fastq_counts.csv', '.fastq_counts.txt'))
            else:
                tab_file_path = str(out_dir / file_path.name.replace('.csv', '.count'))
            
            # Save as tab-delimited
            data.to_csv(tab_file_path, sep="\t", index=False)
            logging.info(f"Created tab-delimited version of count file: {tab_file_path}")
            
        else:
            # This is a tab-delimited file - create CSV version
            tab_file_path = str(file_path)
            # Load the tab-delimited file
            data = pd.read_csv(file_path, sep='\t')
            
            # Create the CSV filename
            if '.counts.txt' in file_path.name.lower():
                csv_file_path = str(out_dir / file_path.name.replace('.counts.txt', '.counts.csv'))
            elif '.fastq_counts.txt' in file_path.name.lower():
                csv_file_path = str(out_dir / file_path.name.replace('.fastq_counts.txt', '.fastq_counts.csv'))
            elif file_path.name.lower().endswith('.count'):
                csv_file_path = str(out_dir / f"{file_path.stem}.counts.csv")
            else:
                csv_file_path = str(out_dir / f"{file_path.stem}.csv")
            
            # Save as CSV
            data.to_csv(csv_file_path, index=False)
            logging.info(f"Created CSV version of count file: {csv_file_path}")
        
        return (csv_file_path, tab_file_path)
        
    except Exception as e:
        logging.error(f"Error ensuring count formats for {count_file_path}: {str(e)}")
        # Return the original file for both formats as a fallback
        return (str(count_file_path), str(count_file_path))


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
    
    # Define patterns for recognizing files
    LIBRARY_PATTERNS = ["*library*.csv", "*library*.txt", "*guide*.csv", "*guide*.txt"]
    
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
    
    # Look for count files directly in the input directory
    cnttbl_path = input_files.get('contrast_file')
    design_path = input_files.get('design_matrix')
    
    # First, check for count files directly in the input directory
    # This is the new recommended location for count files
    logging.info(f"Checking for count files directly in {input_path}")
    
    # Find read count files (RC files)
    rc_files = list(input_path.glob(READ_COUNT_PATTERN))
    
    # Also look for .counts files explicitly
    counts_files = list(input_path.glob("*.counts"))
    
    # Also look for CSV files that match our count patterns
    csv_count_files = []
    for pattern in COUNT_CSV_PATTERN:
        csv_count_files.extend(list(input_path.glob(pattern)))
    
    # Combine all types of count files
    all_count_files = rc_files + counts_files + csv_count_files
    
    if all_count_files and cnttbl_path:
        # Process each count file
        for count_file in all_count_files:
            count_path = str(count_file)
            
            # Ensure the file exists in both CSV and tab-delimited formats
            csv_path, tab_path = ensure_count_formats(count_path)
            
            # Use the tab-delimited version for analysis
            input_files['read_counts'][tab_path] = {
                'contrast_table': cnttbl_path,
                'design_matrix': design_path,
                'csv_path': csv_path  # Store the CSV path for reference
            }
        
        logging.info(f"Found {len(all_count_files)} count files directly in {input_path}")
    
    # Targeted search in standard directories as fallback (for compatibility)
    standard_dirs = ["fastq", "counts", "read_count", "rc"]
    possible_data_dirs = []
    
    # Add standard directories to search paths
    for dir_name in standard_dirs:
        standard_dir = input_path / dir_name
        if standard_dir.exists():
            possible_data_dirs.append(standard_dir)
    
    # Process standard directories if no count files found directly
    if len(input_files['read_counts']) == 0:
        for dir_path in possible_data_dirs:
            if not dir_path.exists():
                continue
                
            logging.info(f"Scanning directory: {dir_path}")
            
            # Find the contrast table (if any) in this directory
            cnttbl_files = list(dir_path.glob(CONTRAST_TABLE_PATTERN))
            cnttbl_path = str(cnttbl_files[0]) if cnttbl_files else input_files.get('contrast_file')
            
            # Find the design matrix (if any) in this directory
            design_files = list(dir_path.glob(DESIGN_MATRIX_PATTERN))
            design_path = str(design_files[0]) if design_files else input_files.get('design_matrix')
            
            # Process potential count files in this directory
            if ("read_count" in str(dir_path).lower() or 
                "_rc" in str(dir_path).lower() or 
                "counts" in str(dir_path).lower() or
                any(f.name.lower().endswith(('.count', '.counts')) for f in dir_path.iterdir() if f.is_file())):
                process_count_directory(dir_path, cnttbl_path, design_path, input_files)
            
            # Check if this is a FASTQ directory
            elif "fastq" in str(dir_path).lower():
                process_fastq_directory(dir_path, cnttbl_path, design_path, input_files)
    
    # If we found very few files, do a full recursive search as fallback
    if len(input_files['read_counts']) == 0 and len(input_files['fastq_dirs']) == 0:
        logging.info("Few files found in standard directories, performing full recursive search...")
        for root, dirs, files in os.walk(input_path):
            root_path = Path(root)
            
            # Skip if no files in directory or already processed
            if not files or str(root_path) in input_files['fastq_dirs']:
                continue
            
            # Find the contrast table (if any) in this directory
            cnttbl_files = [f for f in files if fnmatch_lower(f, CONTRAST_TABLE_PATTERN)]
            cnttbl_path = str(root_path / cnttbl_files[0]) if cnttbl_files else input_files.get('contrast_file')
            
            # Find the design matrix (if any) in this directory
            design_files = [f for f in files if fnmatch_lower(f, DESIGN_MATRIX_PATTERN)]
            design_path = str(root_path / design_files[0]) if design_files else input_files.get('design_matrix')
            
            # Check if this is a read count directory
            if ("read_count" in root.lower() or 
                "_rc" in root.lower() or 
                "counts" in root.lower() or
                any(f.lower().endswith(('.count', '.counts')) for f in files)):
                process_count_directory(root_path, cnttbl_path, design_path, input_files)
            
            # Check if this is a FASTQ directory
            elif "fastq" in root.lower():
                process_fastq_directory(root_path, cnttbl_path, design_path, input_files)
    
    # Count the number of design matrices found
    design_matrices_count = sum(1 for _, info in input_files['read_counts'].items() if info.get('design_matrix') is not None)
    design_matrices_count += sum(1 for _, info in input_files['fastq_dirs'].items() if info.get('design_matrix') is not None)
    
    logging.info(f"Found {len(input_files['read_counts'])} read count files and {len(input_files['fastq_dirs'])} FASTQ directories")
    logging.info(f"Found {design_matrices_count} design matrices for MLE analysis")
    
    return input_files


def process_count_directory(dir_path: Path, cnttbl_path: Optional[str], design_path: Optional[str], input_files: Dict):
    """Process a directory containing count files."""
    # Find read count files (RC files)
    rc_files = list(dir_path.glob(READ_COUNT_PATTERN))
    
    # Also look for .counts files explicitly
    counts_files = list(dir_path.glob("*.counts"))
    
    # Also look for CSV files that match our count patterns
    csv_count_files = []
    for pattern in COUNT_CSV_PATTERN:
        csv_count_files.extend(list(dir_path.glob(pattern)))
    
    # Combine all types of count files
    all_count_files = rc_files + counts_files + csv_count_files
    
    if all_count_files and cnttbl_path:
        # Process each count file
        for count_file in all_count_files:
            count_path = str(count_file)
            
            # Ensure the file exists in both CSV and tab-delimited formats
            csv_path, tab_path = ensure_count_formats(count_path)
            
            # Use the tab-delimited version for analysis
            input_files['read_counts'][tab_path] = {
                'contrast_table': cnttbl_path,
                'design_matrix': design_path,
                'csv_path': csv_path  # Store the CSV path for reference
            }
        
        logging.info(f"Found {len(all_count_files)} count files in {dir_path}")


def process_fastq_directory(dir_path: Path, cnttbl_path: Optional[str], design_path: Optional[str], input_files: Dict):
    """Process a directory containing FASTQ files."""
    # Check for FASTQ files
    fastq_files = []
    for pattern in FASTQ_PATTERNS:
        fastq_files.extend([str(f) for f in dir_path.glob(pattern)])
    
    if fastq_files:
        input_files['fastq_dirs'][str(dir_path)] = {
            'contrast_table': cnttbl_path,
            'design_matrix': design_path,
            'fastq_files': fastq_files
        }
        
        logging.info(f"Found {len(fastq_files)} FASTQ files in {dir_path}")


def ensure_output_dir(output_dir: str) -> None:
    """
    Ensure the output directory exists, creating it if necessary.
    
    Args:
        output_dir: Directory to create
    """
    os.makedirs(output_dir, exist_ok=True)


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


@retry_operation(max_attempts=3, delay=2)
def copy_file(src_path: str, dest_path: str) -> str:
    """
    Copy a file with retry mechanism for file access issues.
    
    Args:
        src_path: Source file path
        dest_path: Destination file path
        
    Returns:
        Path to the copied file
    """
    try:
        logging.info(f"Copying file from {src_path} to {dest_path}")
        shutil.copy2(src_path, dest_path)  # copy2 preserves metadata
        return dest_path
    except Exception as e:
        logging.error(f"Error copying file from {src_path} to {dest_path}: {str(e)}")
        raise


@retry_operation(max_attempts=3, delay=2)
def convert_file_to_tab_delimited(file_path: str, output_dir: Optional[str] = None) -> str:
    """
    Convert a CSV file to tab-delimited format for use with MAGeCK and DrugZ.
    Uses retry_operation decorator to handle temporary file access issues.
    
    Handles both design matrices and contrast tables, including:
    - Cleaning of whitespace in all values
    - Consolidation of duplicate columns in contrast tables
    - Validating required columns for contrast tables
    
    Args:
        file_path: Path to the CSV file
        output_dir: Optional output directory (defaults to same as file_path)
        
    Returns:
        Path to the converted tab-delimited file
    """
    try:
        logging.info(f"Converting {file_path} to tab-delimited format")
        
        # Determine the output path
        input_path = Path(file_path)
        if output_dir:
            output_dir_path = Path(output_dir)
            output_dir_path.mkdir(exist_ok=True, parents=True)
        else:
            output_dir_path = input_path.parent
            
        # Create the output filename (change extension to .txt)
        output_filename = f"{input_path.stem}.txt"
        output_path = output_dir_path / output_filename
        
        # First check if this is a contrast table with duplicate column names
        # Use a more flexible CSV reader for initial inspection
        try:
            with open(file_path, 'r') as f:
                header = f.readline().strip()
            
            # Check for duplicate column names in header
            columns = header.split(',')
            has_duplicates = len(columns) != len(set(columns))
        except Exception as e:
            logging.warning(f"Error reading header to check for duplicates: {str(e)}")
            has_duplicates = False
            
        if has_duplicates:
            logging.info(f"Detected contrast table with duplicate columns")
            
            # Process file with duplicate columns using a custom approach
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            header = lines[0].strip().split(',')
            data_rows = [line.strip().split(',') for line in lines[1:]]
            
            # Find unique column names and their positions
            unique_columns = []
            column_positions = {}
            
            for i, col in enumerate(header):
                col = col.strip()
                if col not in column_positions:
                    column_positions[col] = []
                    unique_columns.append(col)
                column_positions[col].append(i)
            
            # Create a new dataframe with consolidated columns
            new_data = {}
            df_raw = pd.DataFrame(data_rows)
            
            for col in unique_columns:
                positions = column_positions[col]
                if len(positions) > 1:
                    # For duplicate columns, join the values with commas
                    joined_values = []
                    for row in range(len(df_raw)):
                        row_values = []
                        for pos in positions:
                            if pos < len(df_raw.columns):  # Ensure position is valid
                                value = df_raw.iloc[row, pos]
                                if pd.notna(value) and value != '':
                                    # Clean whitespace from the value
                                    cleaned_value = str(value).strip()
                                    if cleaned_value:  # Only add non-empty values
                                        row_values.append(cleaned_value)
                        joined_values.append(','.join(row_values) if row_values else '')
                    new_data[col] = joined_values
                else:
                    # For non-duplicate columns, just copy the values and clean whitespace
                    if positions[0] < len(df_raw.columns):  # Ensure position is valid
                        values = df_raw.iloc[:, positions[0]].values
                        new_data[col] = [str(val).strip() if pd.notna(val) else val for val in values]
            
            # Create the consolidated dataframe
            df = pd.DataFrame(new_data)
            
            # Check if all required columns for contrast tables are present
            missing_columns = [col for col in CONTRAST_REQUIRED_COLUMNS if col not in df.columns]
            if missing_columns:
                raise ValueError(
                    f"Contrast table missing required columns: {', '.join(missing_columns)}. "
                    f"Required columns are: {', '.join(CONTRAST_REQUIRED_COLUMNS)}"
                )
            
            logging.info(f"Consolidated duplicate columns in contrast table")
        else:
            # Try different approaches to read the CSV file
            try:
                # Standard approach first
                df = pd.read_csv(file_path)
            except Exception as e:
                logging.warning(f"Standard CSV reading failed, trying with error_bad_lines=False: {str(e)}")
                try:
                    # Try with on_bad_lines='skip' for newer pandas
                    df = pd.read_csv(file_path, on_bad_lines='skip')
                except TypeError:
                    # Fall back to error_bad_lines for older pandas versions
                    df = pd.read_csv(file_path, error_bad_lines=False)
                
            logging.info(f"Loaded {len(df)} rows from {file_path}")
            
            # Clean whitespace from all string values in the dataframe
            for column in df.columns:
                if df[column].dtype == 'object':  # Only process string/object columns
                    df[column] = df[column].apply(lambda x: x.strip() if isinstance(x, str) else x)
            
            # For contrast tables, clean the comma-separated lists of samples
            is_contrast_table = all(col in df.columns for col in CONTRAST_REQUIRED_COLUMNS)
            if is_contrast_table:
                logging.info(f"File appears to be a contrast table with columns: {', '.join(df.columns)}")
                
                # Clean whitespace from comma-separated values in control and treatment columns
                for col in ['control', 'treatment']:
                    if col in df.columns:
                        df[col] = df[col].apply(lambda x: ','.join([s.strip() for s in str(x).split(',')]) if pd.notna(x) else x)
                
                # Check for required columns
                missing_columns = [col for col in CONTRAST_REQUIRED_COLUMNS if col not in df.columns]
                if missing_columns:
                    raise ValueError(
                        f"Contrast table missing required columns: {', '.join(missing_columns)}. "
                        f"Required columns are: {', '.join(CONTRAST_REQUIRED_COLUMNS)}"
                    )
            else:
                # Assume it's a design matrix
                logging.info(f"File appears to be a design matrix with columns: {', '.join(df.columns)}")
        
        # Save as tab-delimited
        df.to_csv(output_path, sep='\t', index=False)
        logging.info(f"Successfully saved tab-delimited file to {output_path}")
        
        return str(output_path)
        
    except Exception as e:
        logging.error(f"Error converting {file_path} to tab-delimited format: {str(e)}")
        raise 