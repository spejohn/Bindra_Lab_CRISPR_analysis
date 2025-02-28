#!/usr/bin/env python3
"""
Utility script to run DrugZ analysis using Docker.
"""

import os
import subprocess
import argparse
import logging
from pathlib import Path

# Import from other modules
try:
    from analysis_pipeline.core.utils import convert_win_path_to_docker
except ImportError:
    # Try relative import if run as a standalone script
    import sys
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from analysis_pipeline.core.utils import convert_win_path_to_docker

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_drugz_docker(
    count_file: str,
    design_file: str,
    output_dir: str,
    output_prefix: str,
    drugz_args: dict = None
):
    """
    Run DrugZ analysis using Docker.

    Parameters
    ----------
    count_file : str
        Path to the count file.
    design_file : str
        Path to the design file.
    output_dir : str
        Path to the output directory.
    output_prefix : str
        Prefix for output files.
    drugz_args : dict, optional
        Additional arguments to pass to DrugZ, by default None.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get absolute paths
    count_file_abs = os.path.abspath(count_file)
    design_file_abs = os.path.abspath(design_file)
    output_dir_abs = os.path.abspath(output_dir)
    
    # Create directories to mount in Docker
    count_dir = os.path.dirname(count_file_abs)
    design_dir = os.path.dirname(design_file_abs)
    
    # Get filenames for use inside Docker container
    count_filename = os.path.basename(count_file_abs)
    design_filename = os.path.basename(design_file_abs)
    
    # Convert to Docker-compatible paths for Windows
    if os.name == 'nt':  # Windows
        count_dir_docker = convert_win_path_to_docker(count_dir)
        design_dir_docker = convert_win_path_to_docker(design_dir)
        output_dir_docker = convert_win_path_to_docker(output_dir_abs)
        logger.info(f"Running on Windows, converted paths for Docker:")
        logger.info(f"  Count dir: {count_dir} -> {count_dir_docker}")
        logger.info(f"  Design dir: {design_dir} -> {design_dir_docker}")
        logger.info(f"  Output dir: {output_dir_abs} -> {output_dir_docker}")
    else:
        count_dir_docker = count_dir
        design_dir_docker = design_dir
        output_dir_docker = output_dir_abs
    
    # Construct DrugZ command
    drugz_cmd = [
        "python", 
        "/drugz/drugz.py", 
        "-i", f"/count/{count_filename}", 
        "-o", f"/output/{output_prefix}", 
        "-d", f"/design/{design_filename}"
    ]
    
    # Add additional arguments if provided
    if drugz_args:
        for arg, value in drugz_args.items():
            if value is not None:
                if isinstance(value, bool) and value:
                    drugz_cmd.append(f"--{arg}")
                else:
                    drugz_cmd.append(f"--{arg}")
                    drugz_cmd.append(str(value))
    
    # Construct Docker command
    docker_cmd = [
        "docker", "run",
        "--rm",
        "-v", f"{count_dir_docker}:/count",
        "-v", f"{design_dir_docker}:/design",
        "-v", f"{output_dir_docker}:/output",
        "crispr-analysis/drugz:latest"
    ] + drugz_cmd
    
    # Log the command being executed
    logger.info(f"Running DrugZ Docker command: {' '.join(docker_cmd)}")
    logger.info(f"Docker volume mappings:")
    logger.info(f"  Count: {count_dir_docker}:/count")
    logger.info(f"  Design: {design_dir_docker}:/design") 
    logger.info(f"  Output: {output_dir_docker}:/output")
    
    # Run Docker command
    try:
        result = subprocess.run(
            docker_cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.info(f"DrugZ analysis complete: {result.stdout}")
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"DrugZ analysis failed: {e.stderr}")
        return False, e.stderr

def main():
    """
    Main function for running DrugZ analysis using Docker.
    """
    parser = argparse.ArgumentParser(description="Run DrugZ analysis using Docker")
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True, help="Path to count file")
    parser.add_argument("-d", "--design", required=True, help="Path to design file")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output files")
    
    # Optional DrugZ arguments
    parser.add_argument("--minobs", type=int, help="Minimum observations required")
    parser.add_argument("--ctrl", help="Control prefix pattern")
    parser.add_argument("--unpaired", action="store_true", help="Use unpaired analysis")
    parser.add_argument("--normtype", choices=["none", "median", "total"], help="Normalization method")
    parser.add_argument("--build", action="store_true", help="Build index file")
    
    args = parser.parse_args()
    
    # Extract DrugZ-specific arguments
    drugz_args = {
        "minobs": args.minobs,
        "ctrl": args.ctrl,
        "unpaired": args.unpaired,
        "normtype": args.normtype,
        "build": args.build
    }
    
    # Filter out None values
    drugz_args = {k: v for k, v in drugz_args.items() if v is not None}
    
    # Run DrugZ analysis
    success, output = run_drugz_docker(
        args.input,
        args.design,
        args.output,
        args.prefix,
        drugz_args
    )
    
    if success:
        logger.info("DrugZ analysis completed successfully")
    else:
        logger.error("DrugZ analysis failed")
        exit(1)

if __name__ == "__main__":
    main() 