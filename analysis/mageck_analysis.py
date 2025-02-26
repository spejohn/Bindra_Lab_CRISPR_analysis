"""
MAGeCK analysis functions for CRISPR screening.
Implements differential abundance analysis using MAGeCK via Docker.
"""

import os
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any

# Import from other modules
from analysis_pipeline.core.file_handling import ensure_output_dir
from analysis_pipeline.docker.docker_utils import run_docker_container, verify_docker
from analysis_pipeline.core.utils import convert_win_path_to_docker
from analysis_pipeline.core.config import (
    DOCKER_IMAGE,
    DEFAULT_NORM_METHOD,
    DEFAULT_ADJUST_METHOD,
    DEFAULT_FDR_THRESHOLD
)


def mageck_test_analysis(
    count_table: str,
    control_samples: List[str],
    treatment_samples: List[str],
    output_dir: str,
    output_prefix: str,
    norm_method: str = DEFAULT_NORM_METHOD,
    gene_column: str = "Gene",
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD
) -> Tuple[bool, Dict[str, str]]:
    """
    Run MAGeCK test analysis for a single contrast using Docker.
    
    Args:
        count_table: Path to the count table file
        control_samples: List of control sample names
        treatment_samples: List of treatment sample names
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        norm_method: Normalization method (e.g., 'median', 'total')
        gene_column: Column name for gene IDs
        fdr_threshold: FDR threshold for significance
        
    Returns:
        Tuple of (success, output_files_dict)
    """
    # Verify Docker is available
    if not verify_docker():
        error_msg = "Docker verification failed, cannot run MAGeCK test"
        logging.error(error_msg)
        return (False, {"error": error_msg})
        
    # Ensure the output directory exists
    ensure_output_dir(output_dir)
    
    # Prepare Docker volumes
    count_dir = os.path.dirname(os.path.abspath(count_table))
    
    # Create Windows-compatible paths for Docker if needed
    if os.name == 'nt':  # Windows
        volumes = {
            convert_win_path_to_docker(count_dir): {"bind": "/count", "mode": "ro"},
            convert_win_path_to_docker(output_dir): {"bind": "/output", "mode": "rw"}
        }
        logging.info("Running on Windows, converted Docker volume paths: %s", volumes)
    else:
        volumes = {
            count_dir: {"bind": "/count", "mode": "ro"},
            output_dir: {"bind": "/output", "mode": "rw"}
        }
    
    # Convert paths to Docker paths
    docker_count = f"/count/{os.path.basename(count_table)}"
    docker_output_prefix = f"/output/{output_prefix}"
    
    # Join sample names
    control_str = ",".join(control_samples)
    treatment_str = ",".join(treatment_samples)
    
    # Log analysis information
    logging.info(f"Running MAGeCK test analysis:")
    logging.info(f"  Control samples: {control_str}")
    logging.info(f"  Treatment samples: {treatment_str}")
    logging.info(f"  Count table: {count_table}")
    logging.info(f"  Output prefix: {output_prefix}")
    logging.info(f"  Normalization method: {norm_method}")
    
    # Prepare MAGeCK test command
    command = ["mageck", "test"]
    command.extend(["--count-table", docker_count])
    command.extend(["--control-gene", control_str])
    command.extend(["--treatment-gene", treatment_str])
    command.extend(["--output-prefix", docker_output_prefix])
    command.extend(["--norm-method", norm_method])
    command.extend(["--gene-fdr-threshold", str(fdr_threshold)])
    command.extend(["--gene-lfc-method", "median"])
    
    # Enhanced debugging
    logging.info(f"Docker volume mappings: {volumes}")
    logging.info(f"Full MAGeCK test command: {' '.join(command)}")
    
    # Run Docker container
    exit_code, output = run_docker_container(
        command=command,
        volumes=volumes,
        stream_logs=True
    )
    
    if exit_code != 0:
        logging.error(f"MAGeCK test failed for contrast {output_prefix}")
        return (False, {"error": f"MAGeCK test failed: {output}"})
    
    # Define expected output files
    expected_files = {
        "gene_summary": os.path.join(output_dir, f"{output_prefix}.gene_summary.txt"),
        "sgrna_summary": os.path.join(output_dir, f"{output_prefix}.sgrna_summary.txt"),
        "log": os.path.join(output_dir, f"{output_prefix}.log")
    }
    
    # Check if output files exist
    output_files = {}
    missing_files = []
    
    for file_type, file_path in expected_files.items():
        if os.path.exists(file_path):
            output_files[file_type] = file_path
            logging.info(f"MAGeCK output file found: {file_path}")
            
            # Check if the file has content
            if os.path.getsize(file_path) == 0:
                logging.error(f"Output file is empty: {file_path}")
                missing_files.append(file_path)
        else:
            logging.error(f"Expected output file not found: {file_path}")
            missing_files.append(file_path)
    
    if missing_files:
        error_msg = f"MAGeCK test completed but expected output files are missing or empty: {', '.join(missing_files)}"
        logging.error(error_msg)
        output_files["error"] = error_msg
        return (False, output_files)
    
    logging.info(f"MAGeCK test completed successfully for contrast {output_prefix}")
    return (True, output_files)


def process_contrasts(
    contrasts_file: str,
    count_table: str,
    output_dir: str,
    norm_method: str = DEFAULT_NORM_METHOD,
    overwrite: bool = False
) -> Dict[str, Dict[str, str]]:
    """
    Process multiple contrasts defined in a contrasts file.
    
    Args:
        contrasts_file: Path to the contrasts file (CSV)
        count_table: Path to the count table file
        output_dir: Directory for output files
        norm_method: Normalization method
        overwrite: Whether to overwrite existing output files
        
    Returns:
        Dictionary of results for each contrast
    """
    logging.info(f"Processing contrasts from {contrasts_file}")
    
    # Ensure the output directory exists
    ensure_output_dir(output_dir)
    
    # Read contrasts file
    try:
        contrasts_df = pd.read_csv(contrasts_file)
        logging.info(f"Read {len(contrasts_df)} contrasts from {contrasts_file}")
    except Exception as e:
        error_msg = f"Error reading contrasts file {contrasts_file}: {e}"
        logging.error(error_msg)
        return {"error": {"message": error_msg}}
    
    # Validate contrasts file
    required_columns = ["contrast", "control", "treatment"]
    missing_columns = [col for col in required_columns if col not in contrasts_df.columns]
    
    if missing_columns:
        error_msg = f"Contrasts file is missing required columns: {', '.join(missing_columns)}"
        logging.error(error_msg)
        return {"error": {"message": error_msg}}
    
    # Process each contrast
    results = {}
    
    for _, row in contrasts_df.iterrows():
        contrast_name = row["contrast"]
        control_samples = [s.strip() for s in row["control"].split(",")]
        treatment_samples = [s.strip() for s in row["treatment"].split(",")]
        
        logging.info(f"Processing contrast: {contrast_name}")
        logging.info(f"  Control samples: {control_samples}")
        logging.info(f"  Treatment samples: {treatment_samples}")
        
        # Check if output files already exist
        output_prefix = f"{contrast_name}_mageck"
        gene_summary_file = os.path.join(output_dir, f"{output_prefix}.gene_summary.txt")
        
        if os.path.exists(gene_summary_file) and not overwrite:
            logging.info(f"Skipping contrast {contrast_name}, output files already exist")
            results[contrast_name] = {
                "gene_summary": gene_summary_file,
                "sgrna_summary": os.path.join(output_dir, f"{output_prefix}.sgrna_summary.txt"),
                "log": os.path.join(output_dir, f"{output_prefix}.log")
            }
            continue
        
        # Run MAGeCK test
        success, output_files = mageck_test_analysis(
            count_table=count_table,
            control_samples=control_samples,
            treatment_samples=treatment_samples,
            output_dir=output_dir,
            output_prefix=output_prefix,
            norm_method=norm_method
        )
        
        results[contrast_name] = output_files
        
        if not success:
            logging.error(f"MAGeCK test failed for contrast {contrast_name}")
        else:
            logging.info(f"MAGeCK test completed for contrast {contrast_name}")
    
    return results


def run_drugz_analysis(
    count_table: str,
    contrasts_file: str,
    output_dir: str,
    drugz_script_path: Optional[str] = None,
    timeout: int = 900,
    overwrite: bool = False,
    use_docker: bool = True
) -> Dict[str, Dict[str, Any]]:
    """
    Run DrugZ analysis for multiple contrasts.
    
    Args:
        count_table: Path to the count table file
        contrasts_file: Path to the contrasts file (CSV)
        output_dir: Directory for output files
        drugz_script_path: Path to the DrugZ.py script (will be auto-detected if None)
        timeout: Timeout for DrugZ execution in seconds
        overwrite: Whether to overwrite existing output files
        use_docker: Whether to use Docker container for DrugZ (recommended)
        
    Returns:
        Dictionary of results for each contrast
    """
    from analysis_pipeline.core.config import get_drugz_path
    
    logging.info(f"Running DrugZ analysis for contrasts in {contrasts_file}")
    
    # Ensure the output directory exists
    ensure_output_dir(output_dir)

    # Read contrasts file
    try:
        contrasts_df = pd.read_csv(contrasts_file)
        logging.info(f"Read {len(contrasts_df)} contrasts from {contrasts_file}")
    except Exception as e:
        error_msg = f"Error reading contrasts file {contrasts_file}: {e}"
        logging.error(error_msg)
        return {"error": {"message": error_msg}}
    
    # Validate contrasts file
    required_columns = ["contrast", "control", "treatment"]
    missing_columns = [col for col in required_columns if col not in contrasts_df.columns]
    
    if missing_columns:
        error_msg = f"Contrasts file is missing required columns: {', '.join(missing_columns)}"
        logging.error(error_msg)
        return {"error": {"message": error_msg}}
    
    # Process each contrast
    results = {}
    
    for _, row in contrasts_df.iterrows():
        contrast_name = row["contrast"]
        control_samples = [s.strip() for s in row["control"].split(",")]
        treatment_samples = [s.strip() for s in row["treatment"].split(",")]
        
        logging.info(f"Processing DrugZ for contrast: {contrast_name}")
        logging.info(f"  Control samples: {control_samples}")
        logging.info(f"  Treatment samples: {treatment_samples}")
        
        # Check if output files already exist
        output_prefix = f"{contrast_name}_drugz"
        drugz_output_file = os.path.join(output_dir, f"{output_prefix}.txt")
        
        if os.path.exists(drugz_output_file) and not overwrite:
            logging.info(f"Skipping DrugZ for contrast {contrast_name}, output file already exists")
            results[contrast_name] = {"output_file": drugz_output_file}
            continue
        
        # Prepare control and treatment strings
        control_str = ",".join(control_samples)
        treatment_str = ",".join(treatment_samples)

        # Create a temporary design file for this contrast
        design_file = os.path.join(output_dir, f"{contrast_name}_drugz_design.txt")
        try:
            with open(design_file, 'w') as f:
                f.write(f"SAMPLES\tCONDITION\n")
                for sample in control_samples:
                    f.write(f"{sample}\tCTRL\n")
                for sample in treatment_samples:
                    f.write(f"{sample}\tTREAT\n")
            logging.info(f"Created design file for contrast {contrast_name} at {design_file}")
        except Exception as e:
            error_msg = f"Error creating design file for contrast {contrast_name}: {e}"
            logging.error(error_msg)
            results[contrast_name] = {"error": error_msg}
            continue
        
        # Run DrugZ based on use_docker flag
        if use_docker:
            try:
                from analysis_pipeline.analysis.run_drugz_docker import run_drugz_docker
                
                logging.info(f"Running DrugZ in Docker container for contrast {contrast_name}")
                
                # Define arguments for DrugZ
                drugz_args = {
                    "ctrl": "CTRL",
                    "f": "csv",
                    "normtype": "median"
                }
                
                # Run DrugZ in Docker
                success, output = run_drugz_docker(
                    count_file=count_table,
                    design_file=design_file,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    drugz_args=drugz_args
                )
                
                if not success:
                    error_msg = f"DrugZ Docker run failed for contrast {contrast_name}: {output}"
                    logging.error(error_msg)
                    results[contrast_name] = {"error": error_msg}
                    continue
                
                logging.info(f"DrugZ Docker run completed for contrast {contrast_name}")
                
            except Exception as e:
                error_msg = f"Error running DrugZ Docker for contrast {contrast_name}: {e}"
                logging.error(error_msg)
                results[contrast_name] = {"error": error_msg}
                continue
        else:
            # Legacy mode - running without Docker
            # Find DrugZ script if not provided
            if drugz_script_path is None:
                drugz_script_path = get_drugz_path()
                
            if not drugz_script_path or not os.path.exists(drugz_script_path):
                error_msg = f"DrugZ script not found. Please provide the path to DrugZ.py."
                logging.error(error_msg)
                return {"error": {"message": error_msg}}
            
            logging.info(f"Using DrugZ script at {drugz_script_path}")
            
            # Run DrugZ
            try:
                import subprocess
                from subprocess import PIPE
                
                command = [
                    "python", drugz_script_path,
                    "-i", count_table,
                    "-o", drugz_output_file,
                    "-c", control_str,
                    "-t", treatment_str,
                    "-f", "csv"
                ]
                
                logging.info(f"Running DrugZ command: {' '.join(command)}")
                
                # Execute DrugZ with timeout
                process = subprocess.Popen(command, stdout=PIPE, stderr=PIPE)
                
                try:
                    stdout, stderr = process.communicate(timeout=timeout)
                    exit_code = process.returncode
                    
                    stdout_text = stdout.decode() if stdout else ""
                    stderr_text = stderr.decode() if stderr else ""
                    
                    if exit_code != 0:
                        error_msg = f"DrugZ failed with exit code {exit_code} for contrast {contrast_name}: {stderr_text}"
                        logging.error(error_msg)
                        results[contrast_name] = {"error": error_msg}
                        continue
                    
                except subprocess.TimeoutExpired:
                    process.kill()
                    error_msg = f"DrugZ timed out after {timeout} seconds for contrast {contrast_name}"
                    logging.error(error_msg)
                    results[contrast_name] = {"error": error_msg}
                    continue
                    
            except Exception as e:
                error_msg = f"Error running DrugZ for contrast {contrast_name}: {e}"
                logging.error(error_msg)
                results[contrast_name] = {"error": error_msg}
                continue
        
        # Check if output file exists and is not empty
        if not os.path.exists(drugz_output_file):
            error_msg = f"DrugZ output file not found: {drugz_output_file}"
            logging.error(error_msg)
            results[contrast_name] = {"error": error_msg}
            continue
        
        if os.path.getsize(drugz_output_file) == 0:
            error_msg = f"DrugZ output file is empty: {drugz_output_file}"
            logging.error(error_msg)
            results[contrast_name] = {"error": error_msg}
            continue
        
        # Try to read the output file to validate format
        try:
            drugz_df = pd.read_csv(drugz_output_file, sep='\t')
            logging.info(f"DrugZ results file has {len(drugz_df)} rows")
            
            # Convert to CSV for readability
            csv_output_file = os.path.join(output_dir, f"{output_prefix}.csv")
            drugz_df.to_csv(csv_output_file, index=False)
            
            results[contrast_name] = {
                "output_file": drugz_output_file,
                "csv_file": csv_output_file,
                "rows": len(drugz_df)
            }
            
        except Exception as e:
            warning_msg = f"Warning: Error reading DrugZ output file {drugz_output_file}: {e}"
            logging.warning(warning_msg)
            results[contrast_name] = {
                "output_file": drugz_output_file,
                "warning": warning_msg
            }
    
    return results


def find_design_matrix(contrast_dir: str) -> Optional[str]:
    """
    Find a design matrix file in the same directory as the contrasts file.
    
    Args:
        contrast_dir: Directory containing the contrasts file
        
    Returns:
        Path to the design matrix file if found, None otherwise
    """
    logging.info(f"Looking for design matrix file in {contrast_dir}")
    
    # Common names for design matrix files
    design_matrix_patterns = [
        "design_matrix.txt",
        "design.txt",
        "design_matrix.csv",
        "design.csv",
        "design_matrix*.txt",
        "design*.txt"
    ]
    
    # Check for design matrix files
    for pattern in design_matrix_patterns:
        design_files = list(Path(contrast_dir).glob(pattern))
        if design_files:
            design_file = str(design_files[0])
            logging.info(f"Found design matrix file: {design_file}")
            return design_file
    
    logging.info("No design matrix file found")
    return None


def mageck_mle_analysis(
    count_table: str,
    design_matrix: str,
    output_dir: str,
    output_prefix: str,
    norm_method: str = DEFAULT_NORM_METHOD
) -> Tuple[bool, Dict[str, str]]:
    """
    Run MAGeCK MLE analysis using a design matrix via Docker.
    
    Args:
        count_table: Path to the count table file
        design_matrix: Path to the design matrix file
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        norm_method: Normalization method (e.g., 'median', 'total')
        
    Returns:
        Tuple of (success, output_files_dict)
    """
    # Verify Docker is available
    if not verify_docker():
        error_msg = "Docker verification failed, cannot run MAGeCK MLE"
        logging.error(error_msg)
        return (False, {"error": error_msg})
        
    # Ensure the output directory exists
    ensure_output_dir(output_dir)
    
    # Prepare Docker volumes
    count_dir = os.path.dirname(os.path.abspath(count_table))
    design_dir = os.path.dirname(os.path.abspath(design_matrix))
    
    # Create Windows-compatible paths for Docker if needed
    if os.name == 'nt':  # Windows
        volumes = {
            convert_win_path_to_docker(count_dir): {"bind": "/count", "mode": "ro"},
            convert_win_path_to_docker(design_dir): {"bind": "/design", "mode": "ro"},
            convert_win_path_to_docker(output_dir): {"bind": "/output", "mode": "rw"}
        }
        logging.info("Running on Windows, converted Docker volume paths: %s", volumes)
    else:
        volumes = {
            count_dir: {"bind": "/count", "mode": "ro"},
            design_dir: {"bind": "/design", "mode": "ro"},
            output_dir: {"bind": "/output", "mode": "rw"}
        }
    
    # Convert paths to Docker paths
    docker_count = f"/count/{os.path.basename(count_table)}"
    docker_design = f"/design/{os.path.basename(design_matrix)}"
    docker_output_prefix = f"/output/{output_prefix}"
    
    # Log analysis information
    logging.info(f"Running MAGeCK MLE analysis:")
    logging.info(f"  Count table: {count_table}")
    logging.info(f"  Design matrix: {design_matrix}")
    logging.info(f"  Output prefix: {output_prefix}")
    logging.info(f"  Normalization method: {norm_method}")
    
    # Prepare MAGeCK MLE command
    command = ["mageck", "mle"]
    command.extend(["--count-table", docker_count])
    command.extend(["--design-matrix", docker_design])
    command.extend(["--output-prefix", docker_output_prefix])
    command.extend(["--norm-method", norm_method])
    
    # Enhanced debugging
    logging.info(f"Docker volume mappings: {volumes}")
    logging.info(f"Full MAGeCK MLE command: {' '.join(command)}")
    
    # Run Docker container
    exit_code, output = run_docker_container(
        command=command,
        volumes=volumes,
        stream_logs=True
    )
    
    if exit_code != 0:
        logging.error(f"MAGeCK MLE failed for {output_prefix}")
        return (False, {"error": f"MAGeCK MLE failed: {output}"})
    
    # Define expected output files
    expected_files = {
        "gene_summary": os.path.join(output_dir, f"{output_prefix}.gene_summary.txt"),
        "sgrna_summary": os.path.join(output_dir, f"{output_prefix}.sgrna_summary.txt"),
        "log": os.path.join(output_dir, f"{output_prefix}.log")
    }
    
    # Check if output files exist
    output_files = {}
    missing_files = []
    
    for file_type, file_path in expected_files.items():
        if os.path.exists(file_path):
            output_files[file_type] = file_path
            logging.info(f"MAGeCK MLE output file found: {file_path}")
            
            # Check if the file has content
            if os.path.getsize(file_path) == 0:
                logging.error(f"Output file is empty: {file_path}")
                missing_files.append(file_path)
        else:
            logging.error(f"Expected output file not found: {file_path}")
            missing_files.append(file_path)
    
    if missing_files:
        error_msg = f"MAGeCK MLE completed but expected output files are missing or empty: {', '.join(missing_files)}"
        logging.error(error_msg)
        output_files["error"] = error_msg
        return (False, output_files)
    
    logging.info(f"MAGeCK MLE completed successfully for {output_prefix}")
    return (True, output_files)


def process_mle_analysis(
    count_table: str,
    design_matrix: str,
    output_dir: str,
    experiment_name: str,
    norm_method: str = DEFAULT_NORM_METHOD,
    overwrite: bool = False
) -> Dict[str, Dict[str, str]]:
    """
    Process MAGeCK MLE analysis with a design matrix file.
    
    Args:
        count_table: Path to the count table file
        design_matrix: Path to the design matrix file
        output_dir: Directory for output files
        experiment_name: Name of the experiment (used as prefix)
        norm_method: Normalization method
        overwrite: Whether to overwrite existing output files
        
    Returns:
        Dictionary of results
    """
    logging.info(f"Processing MAGeCK MLE analysis with design matrix: {design_matrix}")
    
    # Ensure the output directory exists
    ensure_output_dir(output_dir)
    
    # Check if output files already exist
    output_prefix = f"{experiment_name}_MLE"
    gene_summary_file = os.path.join(output_dir, f"{output_prefix}.gene_summary.txt")
    
    if os.path.exists(gene_summary_file) and not overwrite:
        logging.info(f"Skipping MLE analysis, output files already exist: {gene_summary_file}")
        return {
            "gene_summary": gene_summary_file,
            "sgrna_summary": os.path.join(output_dir, f"{output_prefix}.sgrna_summary.txt"),
            "log": os.path.join(output_dir, f"{output_prefix}.log")
        }
    
    # Run MAGeCK MLE analysis
    success, output_files = mageck_mle_analysis(
        count_table=count_table,
        design_matrix=design_matrix,
        output_dir=output_dir,
        output_prefix=output_prefix,
        norm_method=norm_method
    )
    
    if not success:
        logging.error(f"MAGeCK MLE analysis failed for {experiment_name}")
        return {"error": output_files.get("error", "Unknown error")}
    
    logging.info(f"MAGeCK MLE analysis completed for {experiment_name}")
    return output_files 