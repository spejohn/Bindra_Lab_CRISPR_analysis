"""
Input validation functions for CRISPR analysis pipeline.
"""

import os
import re
import logging
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Union, Tuple, Set

from analysis_pipeline.core.config import (
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