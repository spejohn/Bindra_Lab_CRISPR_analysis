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
    
    # Define expected output files with proper RRA naming convention
    expected_files = {
        "gene_summary": os.path.join(output_dir, f"{output_prefix}_RRA.gene_summary.txt"),
        "sgrna_summary": os.path.join(output_dir, f"{output_prefix}_RRA.sgrna_summary.txt"),
        "log": os.path.join(output_dir, f"{output_prefix}.log")
    }
    
    # Check if MAGeCK output files exist with default names
    default_files = {
        "gene_summary": os.path.join(output_dir, f"{output_prefix}.gene_summary.txt"),
        "sgrna_summary": os.path.join(output_dir, f"{output_prefix}.sgrna_summary.txt"),
        "log": os.path.join(output_dir, f"{output_prefix}.log")
    }
    
    # Rename MAGeCK default files to follow required naming convention
    for file_type in ["gene_summary", "sgrna_summary"]:
        default_file = default_files[file_type]
        expected_file = expected_files[file_type]
        
        if os.path.exists(default_file) and not os.path.exists(expected_file):
            try:
                import shutil
                shutil.move(default_file, expected_file)
                logging.info(f"Renamed {default_file} to {expected_file}")
            except Exception as e:
                logging.error(f"Error renaming {default_file} to {expected_file}: {str(e)}")
    
    # Check if output files exist with expected naming
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
    overwrite: bool = False,
    target_contrast: Optional[str] = None
) -> Dict[str, Dict[str, str]]:
    """
    Process multiple contrasts defined in a contrasts file.
    
    Args:
        contrasts_file: Path to the contrasts file (CSV)
        count_table: Path to the count table file
        output_dir: Directory for output files
        norm_method: Normalization method
        overwrite: Whether to overwrite existing output files
        target_contrast: If provided, only process this specific contrast
        
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
        
        # Skip if this isn't the target contrast
        if target_contrast and contrast_name != target_contrast:
            logging.info(f"Skipping contrast {contrast_name}, not the target contrast")
            continue
            
        control_samples = [s.strip() for s in row["control"].split(",")]
        treatment_samples = [s.strip() for s in row["treatment"].split(",")]
        
        logging.info(f"Processing contrast: {contrast_name}")
        logging.info(f"  Control samples: {control_samples}")
        logging.info(f"  Treatment samples: {treatment_samples}")
        
        # Use the contrast name as the output prefix directly (no need for a subdirectory)
        output_prefix = contrast_name
        gene_summary_file = os.path.join(output_dir, f"{output_prefix}_RRA.gene_summary.txt")
        
        if os.path.exists(gene_summary_file) and not overwrite:
            logging.info(f"Skipping contrast {contrast_name}, output files already exist")
            results[contrast_name] = {
                "gene_summary": gene_summary_file,
                "sgrna_summary": os.path.join(output_dir, f"{output_prefix}_RRA.sgrna_summary.txt"),
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
    use_docker: bool = True,
    target_contrast: Optional[str] = None
) -> Dict[str, Dict[str, Any]]:
    """
    Run DrugZ analysis for all contrasts in the contrasts file.
    
    Args:
        count_table: Path to the count table file
        contrasts_file: Path to the contrasts file
        output_dir: Base directory for output files
        drugz_script_path: Path to the DrugZ script (optional)
        timeout: Timeout for DrugZ analysis in seconds
        overwrite: Whether to overwrite existing output files
        use_docker: Whether to use Docker for running DrugZ
        target_contrast: If provided, only process this specific contrast
        
    Returns:
        Dictionary of results for each contrast
    """
    logging.info(f"Running DrugZ analysis using contrasts from {contrasts_file}")
    
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
        
        # Skip if this isn't the target contrast
        if target_contrast and contrast_name != target_contrast:
            logging.info(f"Skipping contrast {contrast_name}, not the target contrast")
            continue
            
        control_samples = [s.strip() for s in row["control"].split(",")]
        treatment_samples = [s.strip() for s in row["treatment"].split(",")]
        
        logging.info(f"Processing DrugZ for contrast: {contrast_name}")
        logging.info(f"  Control samples: {control_samples}")
        logging.info(f"  Treatment samples: {treatment_samples}")
        
        # Use the contrast name directly for the output file
        output_prefix = contrast_name
        output_file = os.path.join(output_dir, f"{output_prefix}_DrugZ.txt")
        
        if os.path.exists(output_file) and not overwrite:
            logging.info(f"Skipping DrugZ for contrast {contrast_name}, output file already exists")
            results[contrast_name] = {
                "drugz_scores": output_file,
                "log": os.path.join(output_dir, f"{output_prefix}_DrugZ.log")
            }
            continue
        
        # Run DrugZ - implementation will depend on whether using Docker or local installation
        try:
            if use_docker:
                # Paths for Docker
                count_dir = os.path.dirname(os.path.abspath(count_table))
                
                # Create Windows-compatible paths for Docker if needed
                if os.name == 'nt':  # Windows
                    volumes = {
                        convert_win_path_to_docker(count_dir): {"bind": "/count", "mode": "ro"},
                        convert_win_path_to_docker(output_dir): {"bind": "/output", "mode": "rw"}
                    }
                else:
                    volumes = {
                        count_dir: {"bind": "/count", "mode": "ro"},
                        output_dir: {"bind": "/output", "mode": "rw"}
                    }
                
                docker_count = f"/count/{os.path.basename(count_table)}"
                docker_output = f"/output/{output_prefix}_DrugZ.txt"
                
                # Build Docker command
                command = ["python", "/drugz/drugz.py", "-i", docker_count, "-o", docker_output]
                
                # Add control and treatment columns
                for control in control_samples:
                    command.extend(["-c", control])
                for treatment in treatment_samples:
                    command.extend(["-t", treatment])
                
                # Run DrugZ in Docker
                exit_code, output = run_docker_container(
                    image="hartlab/drugz:latest",
                    command=command,
                    volumes=volumes,
                    stream_logs=True,
                    timeout=timeout
                )
                
                if exit_code != 0:
                    logging.error(f"DrugZ analysis failed for contrast {contrast_name}")
                    results[contrast_name] = {"error": f"DrugZ analysis failed: {output}"}
                    continue
                
            else:
                # Run DrugZ locally (simplified example)
                if not drugz_script_path:
                    drugz_script_path = "drugz.py"  # Assume in PATH
                
                # Build command for local execution
                command = [
                    "python", drugz_script_path,
                    "-i", count_table,
                    "-o", output_file
                ]
                
                # Add control and treatment columns
                for control in control_samples:
                    command.extend(["-c", control])
                for treatment in treatment_samples:
                    command.extend(["-t", treatment])
                
                # Run local command (simplified)
                import subprocess
                process = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
                
                if process.returncode != 0:
                    logging.error(f"DrugZ analysis failed for contrast {contrast_name}")
                    results[contrast_name] = {"error": f"DrugZ analysis failed: {process.stderr}"}
                    continue
            
            # Check for output file
            if os.path.exists(output_file):
                logging.info(f"DrugZ analysis completed for contrast {contrast_name}")
                results[contrast_name] = {
                    "drugz_scores": output_file,
                    "log": os.path.join(output_dir, f"{output_prefix}_DrugZ.log")
                }
            else:
                logging.error(f"DrugZ output file not found: {output_file}")
                results[contrast_name] = {"error": "DrugZ output file not found"}
            
        except Exception as e:
            logging.error(f"Error running DrugZ for contrast {contrast_name}: {str(e)}")
            results[contrast_name] = {"error": str(e)}
    
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
    Run MAGeCK MLE analysis using Docker.
    
    Args:
        count_table: Path to the count table file
        design_matrix: Path to the design matrix file
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        norm_method: Normalization method
        
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
        logging.error(f"MAGeCK MLE failed: {output}")
        return (False, {"error": f"MAGeCK MLE failed: {output}"})
    
    # Define expected output files with proper MLE naming
    expected_files = {
        "gene_summary": os.path.join(output_dir, f"{output_prefix}_MLE.gene_summary.txt"),
        "sgrna_summary": os.path.join(output_dir, f"{output_prefix}_MLE.sgrna_summary.txt"),
        "log": os.path.join(output_dir, f"{output_prefix}.log")
    }
    
    # Check if MAGeCK output files exist with default names
    default_files = {
        "gene_summary": os.path.join(output_dir, f"{output_prefix}.gene_summary.txt"),
        "sgrna_summary": os.path.join(output_dir, f"{output_prefix}.sgrna_summary.txt"),
        "log": os.path.join(output_dir, f"{output_prefix}.log")
    }
    
    # Rename MAGeCK default files to follow required naming convention
    for file_type in ["gene_summary", "sgrna_summary"]:
        default_file = default_files[file_type]
        expected_file = expected_files[file_type]
        
        if os.path.exists(default_file) and not os.path.exists(expected_file):
            try:
                import shutil
                shutil.move(default_file, expected_file)
                logging.info(f"Renamed {default_file} to {expected_file}")
            except Exception as e:
                logging.error(f"Error renaming {default_file} to {expected_file}: {str(e)}")
    
    # Check if output files exist with expected naming
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
    
    logging.info(f"MAGeCK MLE completed successfully")
    return (True, output_files)


def process_mle_analysis(
    count_table: str,
    design_matrix: str,
    output_dir: str,
    contrast_name: str,
    norm_method: str = DEFAULT_NORM_METHOD,
    overwrite: bool = False
) -> Dict[str, Dict[str, str]]:
    """
    Process MAGeCK MLE analysis for an experiment.
    
    Args:
        count_table: Path to the count table file
        design_matrix: Path to the design matrix file
        output_dir: Directory for output files
        contrast_name: Name of the contrast
        norm_method: Normalization method
        overwrite: Whether to overwrite existing output files
        
    Returns:
        Dictionary of MLE results
    """
    logging.info(f"Processing MLE analysis for contrast {contrast_name}")
    
    # Ensure the output directory exists
    ensure_output_dir(output_dir)
    
    # Check if output files already exist
    output_prefix = contrast_name
    gene_summary_file = os.path.join(output_dir, f"{output_prefix}_MLE.gene_summary.txt")
    
    if os.path.exists(gene_summary_file) and not overwrite:
        logging.info(f"Skipping MLE analysis, output files already exist")
        return {
            "gene_summary": gene_summary_file,
            "sgrna_summary": os.path.join(output_dir, f"{output_prefix}_MLE.sgrna_summary.txt"),
            "log": os.path.join(output_dir, f"{output_prefix}.log")
        }
    
    # Run MAGeCK MLE
    success, output_files = mageck_mle_analysis(
        count_table=count_table,
        design_matrix=design_matrix,
        output_dir=output_dir,
        output_prefix=output_prefix,
        norm_method=norm_method
    )
    
    if not success:
        logging.error(f"MAGeCK MLE failed for contrast {contrast_name}")
    else:
        logging.info(f"MAGeCK MLE completed for contrast {contrast_name}")
    
    return output_files 