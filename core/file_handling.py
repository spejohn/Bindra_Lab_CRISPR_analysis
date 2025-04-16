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

from core.config import (
    CONTRAST_TABLE_PATTERN,
    READ_COUNT_PATTERN,
    FASTQ_PATTERNS,
    COUNT_CSV_PATTERN,
    DESIGN_MATRIX_PATTERN,
    CONTRAST_REQUIRED_COLUMNS
)
from core.utils import retry_operation

# Define patterns
DESIGN_MATRIX_PATTERN = "*.txt"


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
            with open(file_path, 'r', encoding='utf-8-sig') as f:
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
            with open(file_path, 'r', encoding='utf-8-sig') as f:
                lines = f.readlines()
            
            # Clean header elements thoroughly after splitting
            header = [col.strip() for col in lines[0].strip().split(',')] 
            data_rows = [line.strip().split(',') for line in lines[1:]]
            
            # Find unique column names and their positions
            unique_columns = []
            column_positions = {}
            
            # Use the cleaned header for finding unique columns
            for i, col in enumerate(header):
                # col = col.strip() # Already stripped above
                if col not in column_positions:
                    column_positions[col] = []
                    unique_columns.append(col)
                column_positions[col].append(i)
            
            # Create a new dataframe with consolidated columns
            new_data = {}
            # Create a temporary DataFrame from data_rows to handle varying row lengths safely
            # Ensure columns match the *original* header length for correct indexing
            df_raw = pd.DataFrame(data_rows, columns=header if len(header) == len(data_rows[0]) else None)
            if df_raw.columns is None and len(header) > 0:
                 # If columns couldn't be assigned automatically (ragged array), handle carefully
                 logging.warning("CSV rows have inconsistent lengths. Handling consolidation carefully.")
                 # Reconstruct df_raw column by column to avoid index errors
                 temp_data_for_df = {i: [] for i in range(len(header))}
                 for row in data_rows:
                     for i in range(len(header)):
                         temp_data_for_df[i].append(row[i] if i < len(row) else None)
                 df_raw = pd.DataFrame(temp_data_for_df)
                 df_raw.columns = header # Assign original header names
            
            for col in unique_columns:
                positions = column_positions[col]
                if len(positions) > 1:
                    # For duplicate columns, join the values with commas
                    joined_values = []
                    for row_idx in range(len(df_raw)):
                        row_values = []
                        for pos in positions:
                            if pos < len(df_raw.columns):  # Ensure position is valid
                                value = df_raw.iloc[row_idx, pos]
                                if pd.notna(value) and str(value).strip() != '': # Check for non-empty after stripping
                                    # Clean whitespace from the value
                                    cleaned_value = str(value).strip()
                                    row_values.append(cleaned_value) 
                        joined_values.append(','.join(row_values) if row_values else '') # Use empty string if no valid values
                    new_data[col] = joined_values
                else:
                    # For non-duplicate columns, just copy the values and clean whitespace
                    pos = positions[0]
                    if pos < len(df_raw.columns):  # Ensure position is valid
                        # df_raw might have object dtype, handle potential NAs and types
                        values = df_raw.iloc[:, pos].tolist() # Use tolist() for easier iteration
                        new_data[col] = [str(val).strip() if pd.notna(val) and str(val).strip() != '' else '' for val in values]
                    else:
                         new_data[col] = [''] * len(df_raw) # Add empty column if somehow position is invalid
            
            # Create the consolidated dataframe
            df = pd.DataFrame(new_data)
            # Ensure columns are in the order of unique_columns 
            df = df[unique_columns] 
            
            # Check if all required columns for contrast tables are present *after* consolidation
            # Use the CONTRAST_REQUIRED_COLUMNS list directly
            missing_columns = [col for col in CONTRAST_REQUIRED_COLUMNS if col not in df.columns]
            if missing_columns:
                # Provide more context in the error message
                raise ValueError(
                    f"Contrast table missing required columns after consolidating duplicates: {', '.join(missing_columns)}. "
                    f"Required columns are: {', '.join(CONTRAST_REQUIRED_COLUMNS)}. "
                    f"Columns found after consolidation: {df.columns.tolist()}"
                )
            
            logging.info(f"Consolidated duplicate columns in contrast table. Final columns: {df.columns.tolist()}")
        else:
            # Try different approaches to read the CSV file
            try:
                # Standard approach first, handle BOM with encoding
                df = pd.read_csv(file_path, encoding='utf-8-sig')
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
    
    # Targeted search in standard directories first (more efficient)
    standard_dirs = ["fastq", "counts", "read_count", "rc"]
    possible_data_dirs = [input_path]
    
    # Add standard directories to search paths
    for dir_name in standard_dirs:
        standard_dir = input_path / dir_name
        if standard_dir.exists():
            possible_data_dirs.append(standard_dir)
    
    # Process standard directories first
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
        if "read_count" in str(dir_path).lower() or "_rc" in str(dir_path).lower() or any(f.name.lower().endswith('.count') for f in dir_path.iterdir() if f.is_file()):
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
            if "read_count" in root.lower() or "_rc" in root.lower() or any(f.lower().endswith('.count') for f in files):
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
    
    # Also look for CSV files that match our count patterns
    csv_count_files = []
    for pattern in COUNT_CSV_PATTERN:
        csv_count_files.extend(list(dir_path.glob(pattern)))
    
    # Combine both types of count files
    all_count_files = rc_files + csv_count_files
    
    if all_count_files and cnttbl_path:
        # Process each count file
        for count_file in all_count_files:
            count_path = str(count_file)
            
            # Convert CSV to tab-delimited if needed
            if count_file.name.lower().endswith('.csv'):
                try:
                    count_path = make_count_table(count_path)
                    logging.info(f"Converted CSV count file to tab-delimited: {count_path}")
                except Exception as e:
                    logging.error(f"Error converting CSV count file {count_path}: {str(e)}")
                    continue
            
            input_files['read_counts'][count_path] = {
                'contrast_table': cnttbl_path,
                'design_matrix': design_path
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


def convert_results_to_csv(
    result_file: str, 
    analysis_type: str,
    output_csv_path: Optional[str] = None
) -> str:
    """
    Convert analysis results (RRA, MLE, or DrugZ) to CSV format.
    Selects relevant columns and ensures basic format consistency.
    
    Args:
        result_file: Path to the result file (tab-delimited text output from tool).
        analysis_type: Type of analysis ('RRA', 'MLE', or 'DrugZ').
        output_csv_path: Optional path to save the CSV file. If None, saves next to result_file.
        
    Returns:
        Path to the created CSV file.

    Raises:
        ValueError: If analysis_type is unrecognized or input file is invalid.
        Exception: Re-raises exceptions during file reading or writing.
    """
    logger = logging.getLogger(__name__)
    result_path = Path(result_file)
    analysis_upper = analysis_type.upper()
    logger.info(f"Converting {analysis_upper} results file {result_file} to CSV.")

    if not result_path.is_file():
        raise FileNotFoundError(f"Input result file not found: {result_file}")

    # Determine suffix and expected key columns based on analysis type
    if analysis_upper == 'RRA':
        suffix = '_gMGK'
        # Typical RRA columns: id, ... LFC, p-value, FDR
        # Let's just ensure 'id' (gene ID) is present
        required_col = 'id' 
    elif analysis_upper == 'MLE':
        suffix = '_gMLE'
        # Typical MLE columns: Gene, ... beta, p-value, FDR
        required_col = 'Gene' 
    elif analysis_upper == 'DRUGZ':
        suffix = '_gDZ'
        # Typical DrugZ columns: Gene, ... normZ, pvalue, FDR
        required_col = 'Gene' 
    else:
        raise ValueError(f"Unrecognized analysis_type for CSV conversion: {analysis_type}")
    
    # Determine output path
    if output_csv_path:
        csv_file_path = Path(output_csv_path)
        # Ensure parent directory exists if specified explicitly
        csv_file_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        # Original behavior: save next to input file
        result_dir = result_path.parent
        # Clean up potential existing suffixes in the base name
        base_name = result_path.stem
        for old_suffix in ['.gene_summary', '.drugz', '_RRA', '_MLE', '_DrugZ']:
            if old_suffix in base_name:
                base_name = base_name.replace(old_suffix, '')
        csv_file_path = result_dir / f"{base_name}{suffix}.csv"
    
    try:
        # Read the tab-delimited file, skipping potential comment lines
        df = pd.read_csv(result_file, sep='\t', comment='#')
        
        # Basic validation: Check if the expected key column exists
        if required_col not in df.columns:
             raise ValueError(f"Required column '{required_col}' not found in {result_file}. Found: {df.columns.tolist()}")

        # Select all columns (as per user request - no specific selection/renaming)
        # df_to_save = df[relevant_columns] 
        df_to_save = df
        
        # Save as CSV
        df_to_save.to_csv(csv_file_path, index=False)
        
        logger.info(f"Converted {result_file} to CSV format: {csv_file_path}")
        
        return str(csv_file_path)
    
    except Exception as e:
        logger.error(f"Error converting {result_file} to CSV: {str(e)}")
        # Re-raise the exception to signal failure
        raise Exception(f"Failed to convert {result_file} to CSV.") from e


@retry_operation(max_attempts=3, delay=2)
def copy_file(src_path, dest_path):
    """
    Copy a file from source to destination with error handling and automatic retries.
    
    This function includes robust error handling and uses the retry_operation decorator
    to automatically retry failed copy operations. This is especially useful for
    files that might be temporarily locked or inaccessible due to network issues.
    
    Args:
        src_path: Source file path
        dest_path: Destination file path
        
    Returns:
        Path to the copied file
        
    Raises:
        Exception: If the file could not be copied after multiple attempts
    """
    try:
        # Ensure the destination directory exists
        dest_dir = os.path.dirname(dest_path)
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir, exist_ok=True)
            logging.info(f"Created destination directory: {dest_dir}")
            
        logging.info(f"Copying file from {src_path} to {dest_path}")
        shutil.copy2(src_path, dest_path)  # copy2 preserves metadata
        
        # Verify the file was copied
        if os.path.exists(dest_path):
            logging.info(f"File successfully copied to {dest_path}")
        else:
            logging.error(f"File copy verification failed - destination file does not exist: {dest_path}")
            
        return dest_path
    except Exception as e:
        logging.error(f"Error copying file from {src_path} to {dest_path}: {str(e)}")
        logging.error(f"Source file exists: {os.path.exists(src_path)}")
        logging.error(f"Destination directory exists: {os.path.exists(os.path.dirname(dest_path))}")
        raise


def parse_contrasts(contrasts_txt_path: str) -> List[Dict[str, Union[str, List[str]]]]:
    """
    Parses a tab-delimited contrasts file into a list of dictionaries.

    Expected format (tab-delimited):
    contrast   control         treatment
    exp1_trt   ctrl1,ctrl2     trt1,trt2
    exp2_ko    wt1             ko1,ko2

    Args:
        contrasts_txt_path: Path to the tab-delimited contrast file.

    Returns:
        A list of dictionaries, where each dictionary represents a contrast:
        [{'name': 'exp1_trt', 'control': ['ctrl1','ctrl2'], 'treatment': ['trt1','trt2']},
         {'name': 'exp2_ko', 'control': ['wt1'], 'treatment': ['ko1','ko2']}]

    Raises:
        ValueError: If the file format is incorrect or missing required columns.
        FileNotFoundError: If the file does not exist.
    """
    logger = logging.getLogger(__name__)
    file_path = Path(contrasts_txt_path)

    if not file_path.is_file():
        raise FileNotFoundError(f"Contrast file not found: {contrasts_txt_path}")

    try:
        df = pd.read_csv(file_path, sep='\t')
        logger.info(f"Read contrasts file {file_path} with columns: {df.columns.tolist()}")

        # Check required columns (case sensitive for exact match after conversion)
        required = ['contrast', 'control', 'treatment']
        if not all(col in df.columns for col in required):
            raise ValueError(f"Contrasts file {contrasts_txt_path} must contain columns: {required}. Found: {df.columns.tolist()}")

        parsed_contrasts = []
        for index, row in df.iterrows():
            contrast_name = str(row['contrast']).strip()
            # Split comma-separated samples, strip whitespace, remove empty strings
            control_samples = [s.strip() for s in str(row['control']).split(',') if s.strip()]
            treatment_samples = [s.strip() for s in str(row['treatment']).split(',') if s.strip()]

            if not contrast_name:
                logger.warning(f"Skipping row {index} in {contrasts_txt_path} due to empty contrast name.")
                continue
            if not control_samples:
                logger.warning(f"Skipping contrast '{contrast_name}' in {contrasts_txt_path} due to empty control samples.")
                continue
            if not treatment_samples:
                logger.warning(f"Skipping contrast '{contrast_name}' in {contrasts_txt_path} due to empty treatment samples.")
                continue

            parsed_contrasts.append({
                'name': contrast_name,
                'control': control_samples,
                'treatment': treatment_samples
            })

        if not parsed_contrasts:
             logger.warning(f"No valid contrasts found in {contrasts_txt_path}")

        logger.info(f"Successfully parsed {len(parsed_contrasts)} contrasts from {contrasts_txt_path}")
        return parsed_contrasts

    except Exception as e:
        logger.error(f"Error parsing contrasts file {contrasts_txt_path}: {e}")
        # Re-raise exception after logging
        raise ValueError(f"Error parsing contrasts file {contrasts_txt_path}: {e}") from e 