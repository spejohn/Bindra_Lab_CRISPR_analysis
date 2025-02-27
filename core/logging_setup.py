"""
Logging setup for CRISPR analysis pipeline.
"""

import os
import sys
import logging
import platform
from datetime import datetime
from pathlib import Path


def setup_logging(output_dir: str = None, experiment_name: str = None, log_file: str = None) -> str:
    """
    Configure logging with both file and console output.
    
    Args:
        output_dir: Directory to save log files (optional if log_file is provided)
        experiment_name: Name of the experiment for log file naming (optional if log_file is provided)
        log_file: Direct path to the log file (overrides output_dir and experiment_name)
        
    Returns:
        Path to the log file
    """
    if log_file:
        log_file_path = Path(log_file)
        # Ensure the directory exists
        os.makedirs(log_file_path.parent, exist_ok=True)
    else:
        if not output_dir or not experiment_name:
            raise ValueError("Either log_file or both output_dir and experiment_name must be provided")
            
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%m-%d-%y_%H-%M')
        log_filename = f"{experiment_name}_analysis_{timestamp}.log"
        log_file_path = Path(output_dir) / log_filename

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Log file created at: {log_file_path}")
    return str(log_file_path)


def log_system_info():
    """Log system information to help with debugging."""
    logging.info(f"Python version: {sys.version}")
    logging.info(f"Platform: {platform.platform()}")
    logging.info(f"Python executable: {sys.executable}")
    
    # Log information about installed packages
    try:
        import pandas
        logging.info(f"Pandas version: {pandas.__version__}")
    except ImportError:
        logging.warning("Pandas not installed")
    
    try:
        import numpy
        logging.info(f"NumPy version: {numpy.__version__}")
    except ImportError:
        logging.warning("NumPy not installed")
    
    try:
        import docker
        logging.info(f"Docker SDK version: {docker.__version__}")
    except ImportError:
        logging.warning("Docker SDK not installed")
    
    try:
        import matplotlib
        logging.info(f"Matplotlib version: {matplotlib.__version__}")
    except ImportError:
        logging.warning("Matplotlib not installed")

    # Log Docker version from system
    try:
        import subprocess
        docker_version = subprocess.check_output(["docker", "--version"], stderr=subprocess.STDOUT, text=True)
        logging.info(f"Docker system version: {docker_version.strip()}")
    except (subprocess.SubprocessError, FileNotFoundError):
        logging.warning("Could not determine Docker system version")


def log_input_parameters(parameters: dict):
    """
    Log input parameters for reproducibility.
    
    Args:
        parameters: Dictionary of input parameters
    """
    logging.info("===== Input Parameters =====")
    max_key_length = max(len(str(key)) for key in parameters.keys())
    
    for key, value in parameters.items():
        logging.info(f"{str(key).ljust(max_key_length + 2)}: {value}")
    
    logging.info("=============================") 