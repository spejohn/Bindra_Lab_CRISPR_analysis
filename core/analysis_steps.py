"""
Core analysis step functions, callable independently.
"""

import logging
import os
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union, Any
import pandas as pd

# Assume container_runner is in the same parent directory (core)
try:
    from .container_runner import execute_in_container
except ImportError:
    # Fallback for potential execution context issues
    from container_runner import execute_in_container

# Define default container images (can be overridden via config later)
DEFAULT_FASTQC_IMAGE = "crispr-analysis/fastqc:latest"
DEFAULT_MAGECK_IMAGE = "docker://spejohn/mageck:latest"

logger = logging.getLogger(__name__)


def run_fastqc(
    fastq_path: str,
    output_dir: str,
    threads: int = 1,
    use_apptainer: bool = True,
    container_image: str = DEFAULT_FASTQC_IMAGE,
    timeout: int = 1800 # 30 minutes timeout for FastQC
) -> Tuple[bool, Optional[str]]:
    """
    Runs FastQC on a single FASTQ file using a container.

    Args:
        fastq_path: Path to the input FASTQ file (.fastq.gz or .fq.gz).
        output_dir: Directory where FastQC results should be saved.
        threads: Number of threads for FastQC to use.
        use_apptainer: Whether to prioritize Apptainer for execution.
        container_image: Container image URI (Apptainer SIF or Docker URI).
        timeout: Execution timeout in seconds.

    Returns:
        Tuple (success_boolean, output_report_path_or_error_message).
    """
    logger.info(f"Running FastQC on {fastq_path} using {threads} threads")
    fastq_file = Path(fastq_path)
    host_fastq_dir = fastq_file.parent.resolve()
    host_output_dir = Path(output_dir).resolve()

    # Ensure output directory exists on host
    host_output_dir.mkdir(parents=True, exist_ok=True)

    # Define paths within the container
    container_fastq_dir = "/data/fastq"
    container_output_dir = "/data/output"
    container_fastq_path = f"{container_fastq_dir}/{fastq_file.name}"

    # Define mount points
    mount_map = {
        str(host_fastq_dir): container_fastq_dir,
        str(host_output_dir): container_output_dir
    }

    # Construct FastQC command
    command_list = [
        "fastqc",
        container_fastq_path,
        "-o", container_output_dir,
        "--threads", str(threads)
    ]

    exit_code, output = execute_in_container(
        command_list=command_list,
        container_image=container_image,
        mount_map=mount_map,
        working_dir=container_output_dir, # Run fastqc from output dir
        use_apptainer=use_apptainer,
        timeout=timeout,
        stream_logs=False # Don't stream full FastQC logs by default
    )

    if exit_code == 0:
        # Infer output report name (FastQC usually creates zip and html)
        report_base = fastq_file.name.replace(".fastq.gz", "_fastqc").replace(".fq.gz", "_fastqc")
        html_report = host_output_dir / f"{report_base}.html"
        zip_file = host_output_dir / f"{report_base}.zip"
        
        # Check if expected output files exist
        if html_report.exists() and zip_file.exists():
            logger.info(f"FastQC completed successfully for {fastq_path}. Report: {html_report}")
            return True, str(html_report)
        else:
             logger.warning(f"FastQC exited successfully for {fastq_path}, but expected output files not found ({html_report}, {zip_file}). Log output:\n{output}")
             return False, f"FastQC output files missing. Log:\n{output}"
    else:
        logger.error(f"FastQC failed for {fastq_path} with exit code {exit_code}. Output:\n{output}")
        return False, f"FastQC failed (Exit {exit_code}). Log:\n{output}"

# --- Placeholder for run_mageck_count to be added next ---

def run_mageck_count(
    r1_fastq: str,
    library_path: str,
    output_prefix: str,
    r2_fastq: Optional[str] = None,
    sample_name: Optional[str] = None,
    count_options: Optional[Dict[str, Any]] = None,
    use_apptainer: bool = True,
    container_image: str = DEFAULT_MAGECK_IMAGE,
    timeout: int = 7200 # 2 hours timeout for MAGeCK count
) -> Tuple[bool, Optional[str]]:
    """
    Run MAGeCK count on a single sample (R1, optional R2) using a container.

    Args:
        r1_fastq: Path to the forward reads FASTQ file (.fastq.gz or .fq.gz).
        library_path: Path to the library file (CSV or TXT).
        output_prefix: Prefix for output files (e.g., /path/to/output/samplename).
                       Output files will be samplename.count.txt and samplename.countsummary.txt.
        r2_fastq: Optional path to the reverse reads FASTQ file.
        sample_name: Sample name (defaults to R1 fastq filename stem).
        count_options: Additional options for MAGeCK count (passed as --option value).
        use_apptainer: Whether to prioritize Apptainer for execution.
        container_image: Container image URI (Apptainer SIF or Docker URI).
        timeout: Execution timeout in seconds.

    Returns:
        Tuple (success_boolean, count_file_path_or_error_message).
    """
    logger.info(f"Running MAGeCK count for R1: {r1_fastq}" + (f" R2: {r2_fastq}" if r2_fastq else ""))

    r1_file = Path(r1_fastq)
    r2_file = Path(r2_fastq) if r2_fastq else None
    lib_file = Path(library_path)
    output_pref = Path(output_prefix)

    # Determine sample name if not provided
    if sample_name is None:
        sample_name = r1_file.name.replace(".fastq.gz", "").replace(".fq.gz", "")
    logger.info(f"Using sample name: {sample_name}")

    # Define output paths based on prefix
    host_output_dir = output_pref.parent.resolve()
    output_count_file = output_pref.parent / f"{sample_name}.count.txt"
    output_summary_file = output_pref.parent / f"{sample_name}.countsummary.txt"
    container_output_prefix = f"/data/output/{sample_name}" # Use sample name in prefix

    # Ensure output directory exists on host
    host_output_dir.mkdir(parents=True, exist_ok=True)

    # --- Mount Logic --- 
    # Identify the common directory containing R1 and (potentially) R2 files
    # Assumption: R1 and R2 are in the same directory
    host_fastq_dir = r1_file.parent.resolve()
    host_lib_dir = lib_file.parent.resolve()

    # Define container mount points
    container_fastq_dir = "/data/fastq" 
    container_lib_dir = "/data/library"
    container_output_dir = "/data/output"

    # Define container file paths relative to mount points
    container_r1_path = f"{container_fastq_dir}/{r1_file.name}"
    container_lib_path = f"{container_lib_dir}/{lib_file.name}"
    container_r2_path = f"{container_fastq_dir}/{r2_file.name}" if r2_file else None # Path inside container

    # Define the mount map
    mount_map = {
        str(host_fastq_dir): container_fastq_dir, # Mount the single FASTQ dir
        str(host_lib_dir): container_lib_dir,
        str(host_output_dir): container_output_dir
    }
    # --- End Mount Logic ---

    # Prepare MAGeCK count command
    command_list = ["mageck", "count"]
    command_list.extend(["--fastq", container_r1_path])
    # Use appropriate library file option based on extension
    if lib_file.suffix.lower() == ".csv":
        command_list.extend(["--list-seq", container_lib_path])
    else: # Assume .txt or other format compatible with --library
        command_list.extend(["--library", container_lib_path]) 

    command_list.extend(["--sample-label", sample_name])
    command_list.extend(["--output-prefix", container_output_prefix])

    # Add R2 if provided (path is already defined relative to container_fastq_dir)
    if container_r2_path:
        command_list.extend(["--fastq-2", container_r2_path])

    # Add additional options
    if count_options is None: count_options = {}
    for option, value in count_options.items():
        # Format option keys (remove leading dashes if present)
        clean_option = option.lstrip('-')
        # Append dashes
        formatted_option = f"--{clean_option}"

        if value is None or value is True:
            command_list.append(formatted_option)
        elif value is False:
            pass # Skip false boolean flags
        else:
            command_list.extend([formatted_option, str(value)])

    # Execute in container
    exit_code, output = execute_in_container(
        command_list=command_list,
        container_image=container_image,
        mount_map=mount_map,
        working_dir=container_output_dir,
        use_apptainer=use_apptainer,
        timeout=timeout,
        stream_logs=False # Stream summary later if successful
    )

    if exit_code == 0:
        # Check if expected output count file exists
        if output_count_file.exists():
            logger.info(f"MAGeCK count completed successfully for {sample_name}.")
            logger.info(f"Output count file: {output_count_file}")
            # Optionally log summary content
            if output_summary_file.exists():
                try:
                    summary_content = output_summary_file.read_text()
                    logger.info(f"MAGeCK Count Summary:\n{summary_content}")
                except Exception as e:
                    logger.warning(f"Could not read MAGeCK count summary file {output_summary_file}: {e}")
            else:
                 logger.warning(f"MAGeCK count summary file not found: {output_summary_file}")
            return True, str(output_count_file)
        else:
            logger.error(f"MAGeCK count exited successfully for {sample_name}, but output file {output_count_file} not found. Log output:\n{output}")
            return False, f"MAGeCK count output file missing. Log:\n{output}"
    else:
        logger.error(f"MAGeCK count failed for {sample_name} with exit code {exit_code}. Output:\n{output}")
        return False, f"MAGeCK count failed (Exit {exit_code}). Log:\n{output}"

def aggregate_mageck_counts(
    sample_count_files: List[str],
    output_path: str
) -> Tuple[bool, Optional[str]]:
    """
    Merges individual MAGeCK sample count files into a single aggregated count table.
    Assumes the first two columns are the keys (e.g., sgRNA and Gene) and the third is the count.

    Args:
        sample_count_files: A list of paths to individual sample .count.txt files.
        output_path: The desired path for the aggregated output file.

    Returns:
        Tuple (success_boolean, message_or_output_path).
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Aggregating {len(sample_count_files)} MAGeCK count files into {output_path}")

    if not sample_count_files:
        logger.warning("No sample count files provided for aggregation.")
        # Create an empty file?
        # For now, return True but indicate no files were processed.
        return True, "No sample files provided, no aggregation performed."

    merged_df = None
    key_col_1 = None
    key_col_2 = None

    for i, file_path in enumerate(sample_count_files):
        try:
            df = pd.read_csv(file_path, sep='\t')

            if len(df.columns) < 3:
                 logger.error(f"Count file {file_path} has fewer than 3 columns. Cannot aggregate.")
                 return False, f"File {file_path} has < 3 columns."
            
            current_key_col_1 = df.columns[0]
            current_key_col_2 = df.columns[1]
            current_count_col = df.columns[2]

            # Extract sample name from filename
            sample_name = Path(file_path).stem.replace('.count', '')

            if i == 0:
                # Initialize with the first file
                key_col_1 = current_key_col_1
                key_col_2 = current_key_col_2
                logger.info(f"Using merge keys: '{key_col_1}', '{key_col_2}'")
                # Keep only the key columns and the count column
                merged_df = df[[key_col_1, key_col_2, current_count_col]].copy()
                # Rename the count column to the sample name immediately
                merged_df.rename(columns={current_count_col: sample_name}, inplace=True)
            else:
                # Check if key columns match the first file
                if current_key_col_1 != key_col_1 or current_key_col_2 != key_col_2:
                    logger.error(f"Key columns ('{current_key_col_1}', '{current_key_col_2}') in {file_path} do not match keys ('{key_col_1}', '{key_col_2}') from first file {sample_count_files[0]}. Cannot merge.")
                    return False, "Key column names mismatch between count files."

                # Select keys and the count column
                current_data = df[[key_col_1, key_col_2, current_count_col]].copy()
                # Rename the count column before merging
                current_data.rename(columns={current_count_col: sample_name}, inplace=True)

                # Check if key values/order match before merging (optional but recommended)
                # This ensures rows correspond correctly. Requires keys to be sorted identically.
                # if not merged_df[[key_col_1, key_col_2]].reset_index(drop=True).equals(current_data[[key_col_1, key_col_2]].reset_index(drop=True)):
                #     logger.error(f"Order or content of keys '{key_col_1}', '{key_col_2}' differ between aggregated data and {file_path}. Ensure input files have identical, sorted guides.")
                #     return False, "sgRNA/Gene key content/order mismatch between count files."
                
                # Merge the current sample's count column based on the key columns
                merged_df = pd.merge(merged_df, current_data, on=[key_col_1, key_col_2], how='outer') # Use outer merge to keep all guides

        except Exception as e:
            logger.error(f"Error processing file {file_path}: {e}")
            return False, f"Error reading or processing {file_path}: {e}"

    # Save the final merged dataframe
    try:
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_csv(output_path_obj, sep='\t', index=False, na_rep='NA') # Use NA for missing values
        logger.info(f"Successfully aggregated counts to {output_path}")
        return True, str(output_path_obj)
    except Exception as e:
        logger.error(f"Error saving aggregated count file {output_path}: {e}")
        return False, f"Error saving output file: {e}"

def run_mageck_rra(
    count_path: str,
    contrast_info: Dict[str, Union[str, List[str]]],
    output_prefix: str,
    analysis_options: Optional[Dict[str, Any]] = None,
    use_apptainer: bool = True,
    container_image: str = DEFAULT_MAGECK_IMAGE,
    timeout: int = 3600 # 1 hour timeout for MAGeCK test
) -> Tuple[bool, Optional[str]]:
    """
    Run MAGeCK test (RRA method) for a single contrast using a container.

    Args:
        count_path: Path to the aggregated count file (tab-delimited).
        contrast_info: Dictionary describing the contrast, e.g.,
                       {'name': 'contrast_name', 'control': ['ctrl1'], 'treatment': ['trt1']}.
        output_prefix: Prefix for output files (e.g., /path/to/output/contrast_name).
                       Primary output will be <output_prefix>.gene_summary.txt.
        analysis_options: Additional options for MAGeCK test (passed as --option value).
        use_apptainer: Whether to prioritize Apptainer for execution.
        container_image: Container image URI (Apptainer SIF or Docker URI).
        timeout: Execution timeout in seconds.

    Returns:
        Tuple (success_boolean, gene_summary_path_or_error_message).
    """
    contrast_name = contrast_info.get('name', 'unknown_contrast')
    logger.info(f"Running MAGeCK RRA for contrast: {contrast_name}")

    count_file = Path(count_path)
    output_pref = Path(output_prefix)
    host_output_dir = output_pref.parent.resolve()
    host_count_dir = count_file.parent.resolve()

    # Define expected output file
    output_gene_summary = host_output_dir / f"{output_pref.name}.gene_summary.txt"

    # Ensure output directory exists on host
    host_output_dir.mkdir(parents=True, exist_ok=True)

    # Define container paths
    container_count_dir = "/data/counts"
    container_output_dir = "/data/output"
    container_count_path = f"{container_count_dir}/{count_file.name}"
    container_output_prefix = f"{container_output_dir}/{output_pref.name}"

    # Define mount points
    mount_map = {
        str(host_count_dir): container_count_dir,
        str(host_output_dir): container_output_dir
    }

    # Prepare MAGeCK test command
    command_list = ["mageck", "test"]
    command_list.extend(["-k", container_count_path])

    # Format treatment and control samples
    treatment_samples = contrast_info.get('treatment')
    control_samples = contrast_info.get('control')
    if not treatment_samples or not isinstance(treatment_samples, list):
        err = f"Invalid or missing 'treatment' samples list in contrast_info for {contrast_name}"
        logger.error(err)
        return False, err
    if not control_samples or not isinstance(control_samples, list):
        err = f"Invalid or missing 'control' samples list in contrast_info for {contrast_name}"
        logger.error(err)
        return False, err

    command_list.extend(["-t", ",".join(treatment_samples)])
    command_list.extend(["-c", ",".join(control_samples)])
    command_list.extend(["-n", container_output_prefix])

    # Add additional analysis options
    if analysis_options is None: analysis_options = {}
    # Default norm-method unless specified
    if 'norm-method' not in analysis_options:
        analysis_options['norm-method'] = 'median' 
        
    for option, value in analysis_options.items():
        clean_option = option.lstrip('-')
        formatted_option = f"--{clean_option}"
        if value is None or value is True:
            command_list.append(formatted_option)
        elif value is False:
            pass # Skip false boolean flags
        else:
            command_list.extend([formatted_option, str(value)])

    # Execute in container
    exit_code, output = execute_in_container(
        command_list=command_list,
        container_image=container_image,
        mount_map=mount_map,
        working_dir=container_output_dir,
        use_apptainer=use_apptainer,
        timeout=timeout,
        stream_logs=False # Stream summary later if successful
    )

    if exit_code == 0:
        # Check if expected output gene summary file exists
        if output_gene_summary.exists():
            logger.info(f"MAGeCK RRA completed successfully for contrast {contrast_name}.")
            logger.info(f"Output gene summary: {output_gene_summary}")
            return True, str(output_gene_summary)
        else:
            logger.error(f"MAGeCK RRA exited successfully for {contrast_name}, but output file {output_gene_summary} not found. Log output:\n{output}")
            return False, f"MAGeCK RRA output file missing. Log:\n{output}"
    else:
        logger.error(f"MAGeCK RRA failed for contrast {contrast_name} with exit code {exit_code}. Output:\n{output}")
        return False, f"MAGeCK RRA failed (Exit {exit_code}). Log:\n{output}"

def run_mageck_mle(
    count_path: str,
    design_matrix_path: str,
    output_prefix: str,
    # contrast_info is mainly for naming consistency here, actual comparison derived from design matrix
    contrast_info: Dict[str, Union[str, List[str]]], 
    analysis_options: Optional[Dict[str, Any]] = None,
    use_apptainer: bool = True,
    container_image: str = DEFAULT_MAGECK_IMAGE,
    timeout: int = 7200 # 2 hour timeout for MAGeCK MLE
) -> Tuple[bool, Optional[str]]:
    """
    Run MAGeCK MLE analysis using a design matrix and a container.

    Args:
        count_path: Path to the aggregated count file (tab-delimited).
        design_matrix_path: Path to the design matrix file (tab-delimited).
        output_prefix: Prefix for output files (e.g., /path/to/output/contrast_name_mle).
                       Primary output will be <output_prefix>.gene_summary.txt.
        contrast_info: Dictionary describing the contrast, used mainly for logging/naming.
        analysis_options: Additional options for MAGeCK mle (passed as --option value).
                          Crucially, might include options like --day-label or parameters 
                          defining the model/coefficients based on the design matrix.
        use_apptainer: Whether to prioritize Apptainer for execution.
        container_image: Container image URI (Apptainer SIF or Docker URI).
        timeout: Execution timeout in seconds.

    Returns:
        Tuple (success_boolean, gene_summary_path_or_error_message).
    """
    contrast_name = contrast_info.get('name', 'unknown_contrast')
    logger.info(f"Running MAGeCK MLE for contrast context: {contrast_name}")

    count_file = Path(count_path)
    design_matrix_file = Path(design_matrix_path)
    output_pref = Path(output_prefix)
    host_output_dir = output_pref.parent.resolve()
    host_count_dir = count_file.parent.resolve()
    host_design_dir = design_matrix_file.parent.resolve()

    # Define expected output file
    output_gene_summary = host_output_dir / f"{output_pref.name}.gene_summary.txt"

    # Ensure output directory exists on host
    host_output_dir.mkdir(parents=True, exist_ok=True)

    # Define container paths
    container_count_dir = "/data/counts"
    container_design_dir = "/data/design"
    container_output_dir = "/data/output"
    container_count_path = f"{container_count_dir}/{count_file.name}"
    container_design_path = f"{container_design_dir}/{design_matrix_file.name}"
    container_output_prefix = f"{container_output_dir}/{output_pref.name}"

    # Define mount points
    mount_map = {
        str(host_count_dir): container_count_dir,
        str(host_design_dir): container_design_dir,
        str(host_output_dir): container_output_dir
    }

    # Prepare MAGeCK mle command
    command_list = ["mageck", "mle"]
    command_list.extend(["-k", container_count_path])
    command_list.extend(["-d", container_design_path])
    command_list.extend(["-n", container_output_prefix])

    # Add additional analysis options
    if analysis_options is None: analysis_options = {}
    for option, value in analysis_options.items():
        clean_option = option.lstrip('-')
        formatted_option = f"--{clean_option}"
        if value is None or value is True:
            command_list.append(formatted_option)
        elif value is False:
            pass # Skip false boolean flags
        else:
            command_list.extend([formatted_option, str(value)])

    # Execute in container
    exit_code, output = execute_in_container(
        command_list=command_list,
        container_image=container_image,
        mount_map=mount_map,
        working_dir=container_output_dir,
        use_apptainer=use_apptainer,
        timeout=timeout,
        stream_logs=False # Stream summary later if successful
    )

    if exit_code == 0:
        # Check if expected output gene summary file exists
        if output_gene_summary.exists():
            logger.info(f"MAGeCK MLE completed successfully for contrast context {contrast_name}.")
            logger.info(f"Output gene summary: {output_gene_summary}")
            return True, str(output_gene_summary)
        else:
            logger.error(f"MAGeCK MLE exited successfully for {contrast_name}, but output file {output_gene_summary} not found. Log output:\n{output}")
            return False, f"MAGeCK MLE output file missing. Log:\n{output}"
    else:
        logger.error(f"MAGeCK MLE failed for contrast context {contrast_name} with exit code {exit_code}. Output:\n{output}")
        return False, f"MAGeCK MLE failed (Exit {exit_code}). Log:\n{output}"

# Add default DrugZ image
DEFAULT_DRUGZ_IMAGE = "crispr-analysis/drugz:latest"

def run_drugz(
    count_path: str,
    contrast_info: Dict[str, Union[str, List[str]]],
    output_prefix: str,
    analysis_options: Optional[Dict[str, Any]] = None,
    use_apptainer: bool = True,
    container_image: str = DEFAULT_DRUGZ_IMAGE,
    timeout: int = 3600 # 1 hour timeout for DrugZ
) -> Tuple[bool, Optional[str]]:
    """
    Run DrugZ analysis for a single contrast using a container.

    Args:
        count_path: Path to the aggregated count file (tab-delimited).
        contrast_info: Dictionary describing the contrast, e.g.,
                       {'name': 'contrast_name', 'control': ['ctrl1'], 'treatment': ['trt1']}.
        output_prefix: Prefix for output files (e.g., /path/to/output/contrast_name).
                       Primary output will be <output_prefix>.drugz.txt.
        analysis_options: Additional options for DrugZ script (passed as --option value).
        use_apptainer: Whether to prioritize Apptainer for execution.
        container_image: Container image URI (Apptainer SIF or Docker URI).
        timeout: Execution timeout in seconds.

    Returns:
        Tuple (success_boolean, drugz_output_path_or_error_message).
    """
    contrast_name = contrast_info.get('name', 'unknown_contrast')
    logger.info(f"Running DrugZ for contrast: {contrast_name}")

    count_file = Path(count_path)
    output_pref = Path(output_prefix)
    host_output_dir = output_pref.parent.resolve()
    host_count_dir = count_file.parent.resolve()

    # Define expected output file
    output_drugz_file = host_output_dir / f"{output_pref.name}.drugz.txt"

    # Ensure output directory exists on host
    host_output_dir.mkdir(parents=True, exist_ok=True)

    # Define container paths
    container_count_dir = "/data/counts"
    container_output_dir = "/data/output"
    container_count_path = f"{container_count_dir}/{count_file.name}"
    container_output_path = f"{container_output_dir}/{output_drugz_file.name}"

    # Define mount points
    mount_map = {
        str(host_count_dir): container_count_dir,
        str(host_output_dir): container_output_dir
    }

    # Format treatment and control samples for DrugZ arguments
    treatment_samples = contrast_info.get('treatment')
    control_samples = contrast_info.get('control')
    if not treatment_samples or not isinstance(treatment_samples, list):
        err = f"Invalid or missing 'treatment' samples list in contrast_info for {contrast_name}"
        logger.error(err)
        return False, err
    if not control_samples or not isinstance(control_samples, list):
        err = f"Invalid or missing 'control' samples list in contrast_info for {contrast_name}"
        logger.error(err)
        return False, err

    # Assuming the DrugZ script is at /drugz/drugz.py in the container
    command_list = ["python", "/drugz/drugz.py"]
    command_list.extend(["--input", container_count_path])
    command_list.extend(["--output", container_output_path])
    command_list.extend(["--control-id", ",".join(control_samples)])
    command_list.extend(["--treatment-id", ",".join(treatment_samples)])

    # Add additional analysis options
    if analysis_options is None: analysis_options = {}
    for option, value in analysis_options.items():
        clean_option = option.lstrip('-')
        formatted_option = f"--{clean_option}"
        if value is None or value is True:
            command_list.append(formatted_option)
        elif value is False:
            pass # Skip false boolean flags
        else:
            command_list.extend([formatted_option, str(value)])

    # Execute in container
    exit_code, output = execute_in_container(
        command_list=command_list,
        container_image=container_image,
        mount_map=mount_map,
        working_dir=container_output_dir, # Run from output for simplicity
        use_apptainer=use_apptainer,
        timeout=timeout,
        stream_logs=False
    )

    if exit_code == 0:
        # Check if expected output DrugZ file exists
        if output_drugz_file.exists():
            logger.info(f"DrugZ analysis completed successfully for contrast {contrast_name}.")
            logger.info(f"Output file: {output_drugz_file}")
            return True, str(output_drugz_file)
        else:
            logger.error(f"DrugZ analysis exited successfully for {contrast_name}, but output file {output_drugz_file} not found. Log output:\n{output}")
            return False, f"DrugZ output file missing. Log:\n{output}"
    else:
        logger.error(f"DrugZ analysis failed for contrast {contrast_name} with exit code {exit_code}. Output:\n{output}")
        return False, f"DrugZ analysis failed (Exit {exit_code}). Log:\n{output}"

# --- QC & Reporting Functions --- 

# Attempt to import plotting libraries
try:
    import plotly.express as px
    import plotly.io as pio
    import numpy as np 
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly library not found. Plotting functions will be disabled.")
    # Define placeholder functions if plotly is not available
    def plot_sgRNA_distribution(*args, **kwargs):
        logger.error("Plotly not installed, cannot generate sgRNA distribution plot.")
        return False, "Plotly not installed."
    # Add placeholders for other plot functions here if needed

if PLOTLY_AVAILABLE:
    def plot_sgRNA_distribution(
        count_path: str,
        output_html_path: str,
        log_transform: bool = True,
        pseudocount: int = 1
    ) -> Tuple[bool, Optional[str]]:
        """
        Generates an interactive histogram of sgRNA read count distribution.

        Args:
            count_path: Path to the aggregated count file (tab-delimited).
            output_html_path: Path to save the output HTML file.
            log_transform: Whether to log10 transform counts (+ pseudocount) before plotting.
            pseudocount: Pseudocount to add before log transformation.

        Returns:
            Tuple (success_boolean, output_html_path_or_error_message).
        """
        logger.info(f"Generating sgRNA distribution plot for {count_path}")
        if not PLOTLY_AVAILABLE:
             return False, "Plotly library is required but not installed."

        try:
            df = pd.read_csv(count_path, sep='\t')

            # Identify count columns (all columns except sgRNA and Gene)
            count_cols = [col for col in df.columns if col not in ['sgRNA', 'Gene']]
            if not count_cols:
                raise ValueError("No sample count columns found in the count file.")

            # Calculate total reads per sgRNA across all samples
            # Or plot distribution per sample? Let's start with total per sgRNA
            # Ensure counts are numeric
            df_numeric = df[count_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
            sgRNA_totals = df_numeric.sum(axis=1)
            
            # Create plotting dataframe
            plot_df = pd.DataFrame({'sgRNA': df['sgRNA'], 'counts': sgRNA_totals})
            
            plot_title = "sgRNA Read Count Distribution"
            x_axis_label = "Total Read Counts per sgRNA"
            
            if log_transform:
                plot_df['counts_log10'] = np.log10(plot_df['counts'] + pseudocount)
                x_var = 'counts_log10'
                x_axis_label = f"log10(Counts + {pseudocount})"
                plot_title += f" (log10 + {pseudocount})"
            else:
                x_var = 'counts'

            # Generate histogram using Plotly Express
            fig = px.histogram(
                plot_df,
                x=x_var,
                title=plot_title,
                labels={x_var: x_axis_label},
                marginal="rug", # Add rug plot to show individual points
                # nbins=50 # Optional: adjust number of bins
            )
            fig.update_layout(bargap=0.1)
            
            # Ensure output directory exists
            Path(output_html_path).parent.mkdir(parents=True, exist_ok=True)
            
            # Save plot as HTML
            pio.write_html(fig, output_html_path, auto_open=False)
            
            logger.info(f"Successfully generated sgRNA distribution plot: {output_html_path}")
            return True, str(output_html_path)
            
        except Exception as e:
            logger.exception(f"Error generating sgRNA distribution plot: {e}")
            return False, f"Error generating plot: {e}"

    def plot_gene_distribution(
        count_path: str,
        output_html_path: str,
        log_transform: bool = True,
        pseudocount: int = 1
    ) -> Tuple[bool, Optional[str]]:
        """
        Generates an interactive histogram of gene read count distribution.
        Counts are summed across all sgRNAs targeting the same gene before plotting.

        Args:
            count_path: Path to the aggregated count file (tab-delimited, sgRNA level).
            output_html_path: Path to save the output HTML file.
            log_transform: Whether to log10 transform counts (+ pseudocount) before plotting.
            pseudocount: Pseudocount to add before log transformation.

        Returns:
            Tuple (success_boolean, output_html_path_or_error_message).
        """
        logger.info(f"Generating gene distribution plot for {count_path}")
        # Check is already done at the top level, but double-check for clarity
        if not PLOTLY_AVAILABLE:
            return False, "Plotly library is required but not installed."

        try:
            df = pd.read_csv(count_path, sep='\t')

            # Identify count columns
            count_cols = [col for col in df.columns if col not in ['sgRNA', 'Gene']]
            if not count_cols:
                raise ValueError("No sample count columns found in the count file.")

            # Ensure counts are numeric
            df_numeric = df[count_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
            # Add Gene column back for grouping
            df_numeric['Gene'] = df['Gene']

            # Aggregate counts per gene
            gene_counts = df_numeric.groupby('Gene')[count_cols].sum()
            
            # Calculate total reads per gene across all samples
            gene_totals = gene_counts.sum(axis=1)

            # Create plotting dataframe
            plot_df = pd.DataFrame({'Gene': gene_totals.index, 'counts': gene_totals.values})

            plot_title = "Aggregated Gene Read Count Distribution"
            x_axis_label = "Total Read Counts per Gene"

            if log_transform:
                plot_df['counts_log10'] = np.log10(plot_df['counts'] + pseudocount)
                x_var = 'counts_log10'
                x_axis_label = f"log10(Counts + {pseudocount})"
                plot_title += f" (log10 + {pseudocount})"
            else:
                x_var = 'counts'

            # Generate histogram using Plotly Express
            fig = px.histogram(
                plot_df,
                x=x_var,
                title=plot_title,
                labels={x_var: x_axis_label},
                marginal="rug"
            )
            fig.update_layout(bargap=0.1)

            # Ensure output directory exists
            Path(output_html_path).parent.mkdir(parents=True, exist_ok=True)

            # Save plot as HTML
            pio.write_html(fig, output_html_path, auto_open=False)

            logger.info(f"Successfully generated gene distribution plot: {output_html_path}")
            return True, str(output_html_path)

        except Exception as e:
            logger.exception(f"Error generating gene distribution plot: {e}")
            return False, f"Error generating plot: {e}"

    def plot_gini_index(
        count_path: str,
        output_html_path: str
    ) -> Tuple[bool, Optional[str]]:
        """
        Calculates the Gini index for each sample and generates an interactive Lorenz curve plot.

        Args:
            count_path: Path to the aggregated count file (tab-delimited, sgRNA level).
            output_html_path: Path to save the output HTML file.

        Returns:
            Tuple (success_boolean, output_html_path_or_error_message).
        """
        logger.info(f"Generating Gini index plot for {count_path}")
        if not PLOTLY_AVAILABLE:
            return False, "Plotly library is required but not installed."

        try:
            df = pd.read_csv(count_path, sep='\t')

            # Identify count columns
            count_cols = [col for col in df.columns if col not in ['sgRNA', 'Gene']]
            if not count_cols:
                raise ValueError("No sample count columns found in the count file.")

            # Ensure counts are numeric
            df_numeric = df[count_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

            import plotly.graph_objects as go # Use graph_objects for more control
            fig = go.Figure()
            gini_indices = {}

            # Calculate Gini index and Lorenz curve for each sample
            for sample in count_cols:
                sample_counts = df_numeric[sample].values
                # Remove zero counts as they don't contribute to inequality measure in this context
                sample_counts = sample_counts[sample_counts > 0]
                if len(sample_counts) < 2: # Need at least 2 non-zero counts
                    logger.warning(f"Skipping Gini index for sample '{sample}': Not enough non-zero counts ({len(sample_counts)})." )
                    gini_indices[sample] = np.nan
                    continue
                    
                # Sort counts for Lorenz curve calculation
                sorted_counts = np.sort(sample_counts)
                n = len(sorted_counts)
                cum_counts = np.cumsum(sorted_counts)
                
                # Calculate Lorenz curve points
                lorenz_x = np.linspace(0, 1, n + 1)
                lorenz_y = np.concatenate(([0], cum_counts / cum_counts[-1]))
                
                # Calculate Gini index
                # Gini = 1 - 2 * (Area under Lorenz Curve)
                # Area can be approximated using trapezoidal rule
                area_under_lorenz = np.trapz(lorenz_y, lorenz_x)
                gini = 1 - 2 * area_under_lorenz
                gini_indices[sample] = gini

                # Add Lorenz curve trace to the plot
                fig.add_trace(go.Scatter(
                    x=lorenz_x,
                    y=lorenz_y,
                    mode='lines',
                    name=f'{sample} (Gini: {gini:.3f})',
                    hoverinfo='name+x+y'
                ))

            # Add line of perfect equality
            fig.add_trace(go.Scatter(
                x=[0, 1],
                y=[0, 1],
                mode='lines',
                line=dict(color='grey', dash='dash'),
                name='Perfect Equality'
            ))

            fig.update_layout(
                title='Lorenz Curves and Gini Indices for sgRNA Count Distribution',
                xaxis_title='Cumulative Proportion of sgRNAs',
                yaxis_title='Cumulative Proportion of Reads',
                xaxis=dict(range=[0, 1]),
                yaxis=dict(range=[0, 1]),
                legend_title="Sample (Gini Index)",
                width=800,
                height=700,
                #hovermode='x unified'
            )

            # Ensure output directory exists
            Path(output_html_path).parent.mkdir(parents=True, exist_ok=True)

            # Save plot as HTML
            pio.write_html(fig, output_html_path, auto_open=False)

            logger.info(f"Successfully generated Gini index plot: {output_html_path}")
            # Also return calculated gini indices? For now, just path.
            return True, str(output_html_path)

        except Exception as e:
            logger.exception(f"Error generating Gini index plot: {e}")
            return False, f"Error generating plot: {e}"

    def plot_roc_curve(
        mageck_results_path: str,
        known_controls_path: str,
        output_html_path: str,
        score_col: str = 'rra|score', # Example: Column to use for ranking
        gene_col: str = 'id',
        positive_controls: Optional[List[str]] = None, # Optional list of positive control genes
        negative_controls: Optional[List[str]] = None # Optional list of negative control genes
    ) -> Tuple[bool, Optional[str]]:
        """
        [TODO] Generates an interactive ROC curve based on MAGeCK results and known controls.

        Args:
            mageck_results_path: Path to MAGeCK results file (e.g., gene_summary.txt).
            known_controls_path: Path to file containing lists of known essential (positive) 
                                and non-essential (negative) control genes.
            output_html_path: Path to save the output HTML file.
            score_col: Name of the column in mageck_results_path to use for ranking.
            gene_col: Name of the column containing gene IDs.
            positive_controls: Optional override list of positive control gene names.
            negative_controls: Optional override list of negative control gene names.

        Returns:
            Tuple (success_boolean, output_html_path_or_error_message).
        """
        logger.warning(f"Plotting ROC curve for {mageck_results_path} is not yet implemented.")
        # Placeholder implementation - Needs logic using scikit-learn roc_curve, auc
        # 1. Read MAGeCK results
        # 2. Read known controls file or use provided lists
        # 3. Create binary labels based on known controls
        # 4. Calculate FPR, TPR using roc_curve from sklearn.metrics
        # 5. Calculate AUC
        # 6. Plot using Plotly (go.Scatter)
        # 7. Save HTML
        return False, "ROC curve plotting not implemented."

    # ... [Placeholder/implementation for generate_qc_report goes here] ...


# --- Actual implementations if Plotly is available --- 
if PLOTLY_AVAILABLE:
    # ... [Existing plot_sgRNA_distribution function remains here] ...

    def plot_gene_distribution(
        count_path: str,
        output_html_path: str,
        log_transform: bool = True,
        pseudocount: int = 1
    ) -> Tuple[bool, Optional[str]]:
        """
        Generates an interactive histogram of gene read count distribution.
        Counts are summed across all sgRNAs targeting the same gene before plotting.

        Args:
            count_path: Path to the aggregated count file (tab-delimited, sgRNA level).
            output_html_path: Path to save the output HTML file.
            log_transform: Whether to log10 transform counts (+ pseudocount) before plotting.
            pseudocount: Pseudocount to add before log transformation.

        Returns:
            Tuple (success_boolean, output_html_path_or_error_message).
        """
        logger.info(f"Generating gene distribution plot for {count_path}")
        # Check is already done at the top level, but double-check for clarity
        if not PLOTLY_AVAILABLE:
            return False, "Plotly library is required but not installed."

        try:
            df = pd.read_csv(count_path, sep='\t')

            # Identify count columns
            count_cols = [col for col in df.columns if col not in ['sgRNA', 'Gene']]
            if not count_cols:
                raise ValueError("No sample count columns found in the count file.")

            # Ensure counts are numeric
            df_numeric = df[count_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
            # Add Gene column back for grouping
            df_numeric['Gene'] = df['Gene']

            # Aggregate counts per gene
            gene_counts = df_numeric.groupby('Gene')[count_cols].sum()
            
            # Calculate total reads per gene across all samples
            gene_totals = gene_counts.sum(axis=1)

            # Create plotting dataframe
            plot_df = pd.DataFrame({'Gene': gene_totals.index, 'counts': gene_totals.values})

            plot_title = "Aggregated Gene Read Count Distribution"
            x_axis_label = "Total Read Counts per Gene"

            if log_transform:
                plot_df['counts_log10'] = np.log10(plot_df['counts'] + pseudocount)
                x_var = 'counts_log10'
                x_axis_label = f"log10(Counts + {pseudocount})"
                plot_title += f" (log10 + {pseudocount})"
            else:
                x_var = 'counts'

            # Generate histogram using Plotly Express
            fig = px.histogram(
                plot_df,
                x=x_var,
                title=plot_title,
                labels={x_var: x_axis_label},
                marginal="rug"
            )
            fig.update_layout(bargap=0.1)

            # Ensure output directory exists
            Path(output_html_path).parent.mkdir(parents=True, exist_ok=True)

            # Save plot as HTML
            pio.write_html(fig, output_html_path, auto_open=False)

            logger.info(f"Successfully generated gene distribution plot: {output_html_path}")
            return True, str(output_html_path)

        except Exception as e:
            logger.exception(f"Error generating gene distribution plot: {e}")
            return False, f"Error generating plot: {e}"

    # ... [Placeholders/implementations for plot_gini_index, plot_roc_curve, generate_qc_report go here] ...


# --- Placeholder for plot_gini_index, plot_roc_curve, generate_qc_report --- 

def plot_gini_index(*args, **kwargs):
    logger.error("Plotly not installed, cannot generate Gini index plot.")
    return False, "Plotly not installed."

def plot_roc_curve(*args, **kwargs):
    logger.error("Plotly not installed, cannot generate ROC curve plot.")
    return False, "Plotly not installed."

def generate_qc_report(*args, **kwargs):
    logger.error("Plotly not installed, cannot generate QC report.")
    return False, "Plotly not installed." 