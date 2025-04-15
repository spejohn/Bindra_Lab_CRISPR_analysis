"""
Input validation functions for CRISPR analysis pipeline.
"""

import os
import re
import logging
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Union, Tuple, Set

from core.config import (
    LIBRARY_REQUIRED_COLUMNS,
    CONTRAST_REQUIRED_COLUMNS,
    SAMPLE_SHEET_COLUMNS,
    FASTQ_PATTERNS
)


def validate_library_file(library_file: str, output_dir: str) -> str:
    """
    Validate the library file and potentially create a fixed version.
    
    Args:
        library_file: Path to the library file with guide sequences
        output_dir: Directory to save the fixed library file if needed
        
    Returns:
        Path to the original or fixed library file
    """
    logging.info(f"Validating library file: {library_file}")
    
    try:
        # Try to determine file format
        if library_file.lower().endswith('.csv'):
            df = pd.read_csv(library_file)
            logging.info("Successfully read library file as comma-delimited.")
        else:
            # Try different delimiters
            try:
                df = pd.read_csv(library_file, sep='\t')
                logging.info("Successfully read library file as tab-delimited.")
            except:
                try:
                    df = pd.read_csv(library_file, sep=',')
                    logging.info("Successfully read library file as comma-delimited.")
                except:
                    raise ValueError("Could not determine delimiter for library file.")
        
        logging.info(f"Library file loaded: {len(df)} rows, {len(df.columns)} columns")
        logging.info(f"Columns: {', '.join(df.columns)}")
        
        # Check for required columns
        required_cols = LIBRARY_REQUIRED_COLUMNS
        column_mapping = {}
        needs_fixing = False
        
        # Map column names (case insensitive)
        df_cols_lower = [col.lower() for col in df.columns]
        
        for req_col in required_cols:
            # Check for exact match
            if req_col in df.columns:
                column_mapping[req_col] = req_col
            # Check for case-insensitive match
            elif req_col.lower() in df_cols_lower:
                idx = df_cols_lower.index(req_col.lower())
                orig_col = df.columns[idx]
                column_mapping[req_col] = orig_col
                needs_fixing = True
                logging.info(f"Found column '{orig_col}' matching required '{req_col}'")
            # Check for similar columns
            else:
                similar_cols = []
                for col in df.columns:
                    # Common variations
                    if (req_col.lower() == 'sgrna' and ('guide' in col.lower() or 'grna' in col.lower() or 'sg' in col.lower())):
                        similar_cols.append(col)
                    elif (req_col.lower() == 'gene' and ('target' in col.lower() or 'symbol' in col.lower())):
                        similar_cols.append(col)
                
                if similar_cols:
                    logging.info(f"Required column '{req_col}' not found, but similar columns detected: {similar_cols}")
                    # Use the first similar column
                    column_mapping[req_col] = similar_cols[0]
                    needs_fixing = True
                else:
                    logging.warning(f"Required column '{req_col}' not found and no similar columns detected.")
        
        # Check for missing required columns
        missing_cols = set(required_cols) - set(column_mapping.keys())
        if missing_cols:
            logging.error(f"Missing required columns: {missing_cols}")
            logging.error(f"Available columns: {', '.join(df.columns)}")
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # If we need to fix the library file
        if needs_fixing:
            logging.info("Library file needs fixing. Creating a fixed version...")
            
            # Create a new dataframe with the required columns
            fixed_df = pd.DataFrame()
            
            for req_col, orig_col in column_mapping.items():
                fixed_df[req_col] = df[orig_col]
            
            # Copy any other columns
            for col in df.columns:
                if col not in column_mapping.values():
                    fixed_df[col] = df[col]
            
            # Save the fixed library file
            fixed_path = Path(output_dir) / f"fixed_{Path(library_file).name}"
            if fixed_path.suffix.lower() != '.csv':
                fixed_path = fixed_path.with_suffix('.csv')
                
            fixed_df.to_csv(fixed_path, index=False)
            logging.info(f"Fixed library file saved to: {fixed_path}")
            
            return str(fixed_path)
        else:
            logging.info("Library file appears to be correctly formatted.")
            return library_file
            
    except Exception as e:
        logging.error(f"Error validating library file: {e}")
        raise


def validate_input_tables(
    contrasts: str,
    sample_sheet: Optional[pd.DataFrame] = None,
    counts: Optional[str] = None,
) -> bool:
    """
    Validate input tables for the pipeline.
    
    Args:
        contrasts: Path to tab-delimited file with contrast, treatment, control columns
        sample_sheet: Optional sample sheet DataFrame
        counts: Optional path to count table
        
    Returns:
        True if validation passes, False otherwise
    """
    try:
        # Validate contrasts file
        try:
            contrasts_df = pd.read_csv(Path(contrasts), sep='\t')
        except Exception as e:
            logging.error(f"Error reading contrast table {contrasts}: {str(e)}")
            return False
            
        # Check for required columns
        if set(contrasts_df.columns) != set(CONTRAST_REQUIRED_COLUMNS):
            error_message = (
                f"Error in {contrasts} column names - verify {', '.join(CONTRAST_REQUIRED_COLUMNS)} in table. "
                f"Found columns: {', '.join(contrasts_df.columns)}"
            )
            logging.error(error_message)
            return False
        
        # Validate counts file
        if counts is not None:
            try:
                # First try comma separator
                try:
                    counts_df = pd.read_csv(Path(counts))
                    logging.info(f"Successfully read {counts} with comma separator")
                except:
                    # If that fails, try tab separator
                    try:
                        counts_df = pd.read_csv(Path(counts), sep='\t')
                        logging.info(f"Successfully read {counts} with tab separator")
                    except Exception as e:
                        logging.error(f"Failed to read {counts} with both comma and tab separators: {str(e)}")
                        return False
                
                # Log original columns
                logging.info(f"Original columns in {counts}: {counts_df.columns.tolist()}")
                
                # Convert column names to lowercase and strip whitespace
                columns_lower = {col.lower().strip() for col in counts_df.columns}
                required_columns = {col.lower() for col in LIBRARY_REQUIRED_COLUMNS}
                
                # Log processed columns
                logging.info(f"Processed lowercase columns: {sorted(list(columns_lower))}")
                
                # Check for similar column names that might be variants
                for col in counts_df.columns:
                    if 'guide' in col.lower() or 'sgrna' in col.lower():
                        logging.info(f"Found potential guide column: {col}")
                    if 'gene' in col.lower():
                        logging.info(f"Found potential gene column: {col}")
                
                if not required_columns.issubset(columns_lower):
                    missing_cols = required_columns - columns_lower
                    error_message = (
                        f"{counts} is missing required columns: {missing_cols}.\n"
                        f"Found columns (original): {counts_df.columns.tolist()}\n"
                        f"Found columns (lowercase): {sorted(list(columns_lower))}\n"
                        f"First few rows:\n{counts_df.head().to_string()}"
                    )
                    logging.error(error_message)
                    return False
                    
                logging.info(f"Successfully validated count table {counts}")
                
            except Exception as e:
                error_message = f"Error reading count table {counts}: {str(e)}"
                logging.error(error_message)
                return False

        # Validate sample sheet
        if sample_sheet is not None:
            try:
                required_columns = set(SAMPLE_SHEET_COLUMNS)
                if set(sample_sheet.columns) != required_columns:
                    error_message = (
                        f"Sample sheet columns are expected to be exactly: {', '.join(required_columns)}. "
                        f"Found columns: {', '.join(sample_sheet.columns)}"
                    )
                    logging.error(error_message)
                    return False
                
                # Verify that all files in the sample sheet exist
                for _, row in sample_sheet.iterrows():
                    fastq_paths = row['fastq_path'].split(',')
                    for path in fastq_paths:
                        if not Path(path).exists():
                            logging.warning(f"FASTQ file does not exist: {path}")
                
                logging.info(f"Successfully validated sample sheet with {len(sample_sheet)} samples")
            except Exception as e:
                error_message = f"Error validating sample sheet: {str(e)}"
                logging.error(error_message)
                return False
            
        return True
        
    except Exception as e:
        logging.error(f"Error in validate_input_tables: {str(e)}")
        return False


def check_existing_files(input_dir: str, output_dir: str) -> Set[str]:
    """
    Check for existing analysis files to avoid reprocessing.
    
    Args:
        input_dir: Input directory path
        output_dir: Output directory path
        
    Returns:
        Set of file paths that have already been analyzed
    """
    input_path = Path(input_dir)
    analyzed_files = set()  # Track specific files that have been analyzed

    try:
        # Iterate through input_dir for files ending in "_RC" or ".fastq"
        for root, dirs, files in os.walk(input_path):
            for file in files:
                # Skip xlsx files
                if file.endswith('.xlsx'):
                    continue
                    
                if "_rc" in file.lower() or any(file.lower().endswith(ext) for ext in ['.fq', '.fastq', '.fastq.gz', '.fq.gz']):
                    try:
                        # Create the corresponding output path
                        dir_str = re.sub(f'^{re.escape(str(input_path))}', '', root)
                        existing_output = Path(output_dir) / Path(dir_str.lstrip('\\'))
                        
                        # Check for MAGeCK and DrugZ output files
                        mgk_files = list(existing_output.glob("*_gMGK.csv"))
                        dz_files = list(existing_output.glob("*_gDZ.csv"))

                        if mgk_files or dz_files:
                            logging.info(f"Analysis for {file} exists in {existing_output}")
                            analyzed_files.add(str(Path(root) / file))
                    except Exception as e:
                        logging.warning(f"Error checking output for {file}: {str(e)}")
                
    except Exception as e:
        logging.error(f"Error checking existing files: {str(e)}")
    
    logging.info(f"Found {len(analyzed_files)} previously analyzed files")
    return analyzed_files


def validate_experiment_structure(experiment_dir: str) -> Dict[str, Union[str, bool, List[str]]]:
    """
    Checks the existence and basic structure of required files within an experiment directory.

    Args:
        experiment_dir: Path to the experiment directory.

    Returns:
        A dictionary summarizing findings:
        {'status': 'valid',
         'data_type': 'fastq' or 'rc',
         'contrasts_path': path_to_contrasts.csv,
         'library_path': path_to_library.csv or None,
         'rc_path': path_to_rc_file or None,
         'design_matrix_path': path_to_design_matrix.csv or None,
         'fastq_dir': path_to_fastq_dir or None,
         'fastq_files': list_of_fastq_files or None,
         'mle_possible': True or False}

    Raises:
        ValueError: If the directory structure is invalid or missing required files.
    """
    logger = logging.getLogger(__name__)
    exp_path = Path(experiment_dir)
    logger.info(f"Validating structure of experiment directory: {exp_path}")

    if not exp_path.is_dir():
        raise ValueError(f"Experiment directory not found: {experiment_dir}")

    # --- Check for required contrasts file --- 
    # Allow for flexible naming patterns (*contrasts.csv/txt, *_cnttbl.csv/txt)
    contrast_patterns = ["*contrasts.csv", "*_cnttbl.csv", "*contrasts.txt", "*_cnttbl.txt"]
    contrasts_files = []
    for pattern in contrast_patterns:
        contrasts_files.extend(list(exp_path.glob(pattern)))
    
    # Filter out directories if any match the pattern (glob can sometimes include dirs)
    contrasts_files = [f for f in contrasts_files if f.is_file()]

    if not contrasts_files:
        raise ValueError(f"Missing required contrasts file (matching patterns: {contrast_patterns}) in {experiment_dir}")
    if len(contrasts_files) > 1:
        # Sort for deterministic selection if multiple found
        contrasts_files.sort()
        logger.warning(f"Multiple contrast files found matching patterns in {experiment_dir}. Using the first one: {contrasts_files[0]}")
    contrasts_path = str(contrasts_files[0])
    logger.info(f"Found contrasts file: {contrasts_path}")

    # --- Check for data type (FASTQ or RC) --- 
    fastq_dir = exp_path / "fastq"
    # Look for read count files (e.g., *_rc.csv, experiment_counts.csv)
    rc_files = list(exp_path.glob("*_rc.csv")) + list(exp_path.glob("experiment_counts.csv"))

    has_fastq = fastq_dir.is_dir()
    has_rc = bool(rc_files)

    data_type = None
    library_path = None
    rc_path = None
    fastq_dir_path = None
    fastq_files_list = None

    if has_fastq:
        data_type = "fastq"
        fastq_dir_path = str(fastq_dir)
        fastq_files_list = [str(f) for f in fastq_dir.glob("*.fastq.gz")] + \
                           [str(f) for f in fastq_dir.glob("*.fq.gz")]
        if not fastq_files_list:
             raise ValueError(f"fastq directory found in {experiment_dir}, but it contains no .fastq.gz or .fq.gz files.")
        logger.info(f"Found fastq data directory: {fastq_dir_path} with {len(fastq_files_list)} fastq.gz/fq.gz files.")

        # Library file is required for FASTQ
        library_files = list(exp_path.glob("library.csv"))
        if not library_files:
            raise ValueError(f"Missing required library file (library.csv) for FASTQ processing in {experiment_dir}")
        if len(library_files) > 1:
             logger.warning(f"Multiple library files found in {experiment_dir}, using {library_files[0]}")
        library_path = str(library_files[0])
        logger.info(f"Found library file: {library_path}")

        if has_rc:
            logger.warning(f"Both fastq directory and read count file(s) ({[str(f) for f in rc_files]}) found in {experiment_dir}. Prioritizing fastq for analysis type determination.")
            # Decide on priority or specific handling TBD - for now, just note it.

    elif has_rc:
        data_type = "rc"
        if len(rc_files) > 1:
            logger.warning(f"Multiple read count files found in {experiment_dir}, using {rc_files[0]}")
        rc_path = str(rc_files[0])
        logger.info(f"Found read count file: {rc_path}")
    else:
        raise ValueError(f"Missing data files in {experiment_dir}. Expected either a 'fastq/' subdirectory or a '*_rc.csv'/'experiment_counts.csv' file.")

    # --- Check for optional design matrix --- 
    # Allow for .csv or .txt initially
    design_files = list(exp_path.glob("design_matrix.csv")) + list(exp_path.glob("design_matrix.txt"))
    design_matrix_path = None
    mle_possible = False
    if design_files:
        if len(design_files) > 1:
             logger.warning(f"Multiple design matrix files found in {experiment_dir}, using {design_files[0]}")
        design_matrix_path = str(design_files[0])
        mle_possible = True
        logger.info(f"Found design matrix file: {design_matrix_path} (MLE possible)")
    else:
        logger.info(f"No design matrix file found in {experiment_dir} (MLE not possible).")

    # --- Return summary --- 
    validation_summary = {
        'status': 'valid',
        'data_type': data_type,
        'contrasts_path': contrasts_path,
        'library_path': library_path,
        'rc_path': rc_path,
        'design_matrix_path': design_matrix_path,
        'fastq_dir': fastq_dir_path,
        'fastq_files': fastq_files_list,
        'mle_possible': mle_possible
    }
    logger.info(f"Validation summary for {experiment_dir}: {validation_summary}")
    return validation_summary 