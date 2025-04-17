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
    Validates the structure of a single experiment directory.
    Checks for required files (contrasts, library) and identifies data type.

    Args:
        experiment_dir: Path to the experiment directory.

    Returns:
        A dictionary containing validation status, paths, and detected data type.
        Example: {'status': 'valid', 'data_type': 'fastq', 'contrasts_path': '...', ...}
               {'status': 'failed', 'error': 'Reason for failure'}
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Validating structure of experiment directory: {experiment_dir}")
    exp_path = Path(experiment_dir)
    info: Dict[str, Union[str, bool, List[str], None]] = {
        "status": "valid",
        "data_type": None,
        "contrasts_path": None,
        "library_path": None,
        "rc_path": None,
        "design_matrix_path": None,
        "fastq_dir": None,
        "fastq_files": [],
        "mle_possible": False,
        "error": None
    }

    if not exp_path.is_dir():
        err_msg = f"Experiment directory not found: {experiment_dir}"
        logger.error(err_msg)
        info["status"] = "failed"
        info["error"] = err_msg
        return info

    # --- Find Required Files --- 

    # 1. Find Contrasts File (mandatory)
    contrast_files = list(exp_path.glob("*contrast*.csv"))
    contrast_files.extend(list(exp_path.glob("*contrast*.txt")))
    if not contrast_files:
        err_msg = f"Mandatory contrasts file (*contrast*.csv or *contrast*.txt) not found in {experiment_dir}."
        logger.error(err_msg)
        info["status"] = "failed"
        info["error"] = err_msg
        return info
    elif len(contrast_files) > 1:
        logger.warning(f"Multiple contrast files found in {experiment_dir}. Using the first one: {contrast_files[0]}")
    info["contrasts_path"] = str(contrast_files[0])
    logger.info(f"Found contrasts file: {info['contrasts_path']}")

    # 2. Find Library File (mandatory)
    library_patterns = ["*library*.csv", "*library*.txt", "*guide*.csv", "*guide*.txt"]
    library_files = []
    for pattern in library_patterns:
        library_files.extend(list(exp_path.glob(pattern)))
    
    if not library_files:
        err_msg = f"Mandatory library/guide file ({'/'.join(library_patterns)}) not found in {experiment_dir}."
        logger.error(err_msg)
        info["status"] = "failed"
        info["error"] = err_msg
        return info
    elif len(library_files) > 1:
         logger.warning(f"Multiple library/guide files found in {experiment_dir}. Using the first one: {library_files[0]}")
    info["library_path"] = str(library_files[0])
    logger.info(f"Found library file: {info['library_path']}")
    
    # *** Add Library Header Order Check ***
    try:
        # Read only the header using pandas to auto-detect separator
        lib_header_df = pd.read_csv(
            info["library_path"], 
            sep=None, # Auto-detect separator (comma or tab)
            engine='python', # Needed for sep=None
            nrows=0, # Read only header
            encoding='utf-8-sig' # Handle potential BOM
        )
        lib_columns = [col.lower().strip() for col in lib_header_df.columns]
        
        # Use LIBRARY_REQUIRED_COLUMNS (lowercase) for order check
        expected_order_lower = [col.lower() for col in LIBRARY_REQUIRED_COLUMNS]
        
        if len(lib_columns) < 3:
             # Use expected_order_lower in the error message
            raise ValueError(f"Library file has fewer than 3 columns ({len(lib_columns)} found). Expected at least: {expected_order_lower}")
            
        # Check the first three columns in order (case-insensitive)
        first_three_cols = lib_columns[:3]
        # Compare against expected_order_lower
        if first_three_cols != expected_order_lower:
             raise ValueError(
                f"Incorrect library file column order. Expected first three columns: "
                # Use expected_order_lower in the error message
                f"{expected_order_lower}. Found: {first_three_cols}. "
                f"Full columns found: {lib_columns}"
            )
        logger.info(f"Library file header order validated successfully: {first_three_cols}...")
            
    except Exception as e:
        err_msg = f"Error validating library file header ({info['library_path']}): {e}"
        logger.error(err_msg)
        info["status"] = "failed"
        info["error"] = err_msg
        return info
    # *** End Library Header Order Check ***

    # --- Identify Data Type (FASTQ or Read Counts) --- 
    # Prefer FASTQ if found
    fastq_dir_path = exp_path / "fastq"
    has_fastq_files = False
    if fastq_dir_path.is_dir():
        for pattern in FASTQ_PATTERNS:
            found_files = list(fastq_dir_path.glob(pattern))
            if found_files:
                info["fastq_files"].extend([str(f) for f in found_files])
                has_fastq_files = True
        if has_fastq_files:
            info["data_type"] = "fastq"
            info["fastq_dir"] = str(fastq_dir_path)
            logger.info(f"Found fastq data directory: {info['fastq_dir']} with {len(info['fastq_files'])} fastq.gz/fq.gz files.")

    # If no FASTQ, check for read count file
    if not has_fastq_files:
        rc_patterns = ["*read_count*.csv", "*read_count*.txt", "*counts*.csv", "*counts*.txt"]
        rc_files = []
        for pattern in rc_patterns:
             # Search directly in experiment_dir first
            rc_files.extend(list(exp_path.glob(pattern)))
             # Then check common subdirs like 'counts' or 'rc'
            for sub in ["counts", "rc", "read_counts"]:
                rc_files.extend(list((exp_path / sub).glob(pattern)))
        
        # Filter out library/contrast files that might match count patterns
        rc_files = [
            f for f in rc_files 
            if f.name != Path(info["contrasts_path"]).name and 
               f.name != Path(info["library_path"]).name
        ]
        
        if rc_files:
            info["data_type"] = "rc" # read counts
            info["rc_path"] = str(rc_files[0]) # Use the first one found
            if len(rc_files) > 1:
                logger.warning(f"Multiple read count files found in {experiment_dir} or subdirs. Using the first one: {info['rc_path']}")
            logger.info(f"Found read count data file: {info['rc_path']}")
        else:
            # If neither FASTQ nor RC file found
            err_msg = f"No FASTQ data directory or read count file found in {experiment_dir}. Cannot determine data type."
            logger.error(err_msg)
            info["status"] = "failed"
            info["error"] = err_msg
            return info

    # --- Optional: Find Design Matrix for MLE --- 
    design_matrix_patterns = ["*design_matrix*.csv", "*design_matrix*.txt", "*design*.csv", "*design*.txt"]
    design_files = []
    for pattern in design_matrix_patterns:
        design_files.extend(list(exp_path.glob(pattern)))
    
    # Filter out contrast/library/rc files that might match design patterns
    potential_conflicts = {Path(p).name for p in [info["contrasts_path"], info["library_path"], info["rc_path"]] if p}
    design_files = [f for f in design_files if f.name not in potential_conflicts]

    if design_files:
        info["design_matrix_path"] = str(design_files[0])
        info["mle_possible"] = True
        if len(design_files) > 1:
            logger.warning(f"Multiple design matrix files found in {experiment_dir}. Using the first one: {info['design_matrix_path']}")
        logger.info(f"Found design matrix file: {info['design_matrix_path']} (MLE possible)")
    else:
        logger.info(f"No design matrix file found in {experiment_dir} (MLE not possible).")
        info["mle_possible"] = False

    # Final status check
    if info["status"] == "valid":
        logger.info(f"Validation summary for {experiment_dir}: {info}")
    else:
        # Log the error again for clarity before returning
        logger.error(f"Validation failed for {experiment_dir}: {info['error']}")
        
    return info


# Example usage (if run directly)
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # Create dummy structure for testing
    test_dir = Path("temp_validation_test")
    test_dir.mkdir(exist_ok=True)
    (test_dir / "exp1").mkdir(exist_ok=True)
    (test_dir / "exp1" / "my_contrasts.csv").touch()
    (test_dir / "exp1" / "my_library.csv").write_text("sgrna,sequence,gene\nAG1,AAA,GENE1\nAG2,CCC,GENE2")
    (test_dir / "exp1" / "fastq").mkdir(exist_ok=True)
    (test_dir / "exp1" / "fastq" / "sample1_R1.fastq.gz").touch()
    (test_dir / "exp1" / "fastq" / "sample1_R2.fastq.gz").touch()
    (test_dir / "exp1" / "exp1_design.txt").touch()

    (test_dir / "exp2").mkdir(exist_ok=True)
    (test_dir / "exp2" / "my_contrasts_exp2.txt").touch()
    (test_dir / "exp2" / "my_guide_library.txt").write_text("sgrna\tgene\tsequence\nAG1\tGENE1\tAAA\nAG2\tGENE2\tCCC") # Incorrect order
    (test_dir / "exp2" / "counts").mkdir(exist_ok=True)
    (test_dir / "exp2" / "counts" / "exp2_counts.csv").touch()
    
    (test_dir / "exp3_nolib").mkdir(exist_ok=True)
    (test_dir / "exp3_nolib" / "exp3_contrasts.csv").touch()
    (test_dir / "exp3_nolib" / "fastq").mkdir(exist_ok=True)
    (test_dir / "exp3_nolib" / "fastq" / "s1.fq.gz").touch()

    print("--- Testing exp1 (FASTQ, Design Matrix, Correct Lib Header) ---")
    info1 = validate_experiment_structure(str(test_dir / "exp1"))
    print(info1)
    assert info1["status"] == "valid"
    assert info1["data_type"] == "fastq"
    assert info1["mle_possible"] == True

    print("\n--- Testing exp2 (RC, Incorrect Lib Header) ---")
    info2 = validate_experiment_structure(str(test_dir / "exp2"))
    print(info2)
    assert info2["status"] == "failed"
    assert "Incorrect library file column order" in info2["error"]

    print("\n--- Testing exp3 (No Library) ---")
    info3 = validate_experiment_structure(str(test_dir / "exp3_nolib"))
    print(info3)
    assert info3["status"] == "failed"
    assert "Mandatory library/guide file" in info3["error"]

    # Clean up dummy files
    import shutil
    shutil.rmtree(test_dir)
    print("\nCleanup complete.") 