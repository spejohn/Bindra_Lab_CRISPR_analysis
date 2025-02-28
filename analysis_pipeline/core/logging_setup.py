"""
Logging setup for CRISPR analysis pipeline.
"""

import os
import sys
import logging
import time
from typing import Dict, Any, Optional
import platform
from datetime import datetime
from pathlib import Path

class ProgressReporter:
    """
    Progress reporting tool for tracking pipeline progress.
    Provides both logging and console output.
    """
    
    def __init__(self, total_steps: int, experiment_name: str, tty_output: bool = True):
        """
        Initialize the progress reporter.
        
        Args:
            total_steps: Total number of steps in the pipeline
            experiment_name: Name of the experiment
            tty_output: Whether to output to console (TTY)
        """
        self.total_steps = total_steps
        self.current_step = 0
        self.experiment_name = experiment_name
        self.start_time = time.time()
        self.last_update_time = self.start_time
        self.tty_output = tty_output
        self.completed_steps = []
        
        # Initial message
        self._log_status(f"Starting pipeline for {experiment_name}")
    
    def update(self, step_name: str, step_description: Optional[str] = None):
        """
        Update progress.
        
        Args:
            step_name: Name of the completed step
            step_description: Optional description of the step
        """
        self.current_step += 1
        self.completed_steps.append(step_name)
        current_time = time.time()
        
        elapsed = current_time - self.start_time
        step_time = current_time - self.last_update_time
        percent = (self.current_step / self.total_steps) * 100
        
        # Build status message
        status = f"Progress: {percent:.1f}% | Step {self.current_step}/{self.total_steps}: {step_name}"
        if step_description:
            status += f" - {step_description}"
        status += f" (step time: {step_time:.1f}s, total elapsed: {elapsed:.1f}s)"
        
        # Log status
        self._log_status(status)
        
        # Update last time
        self.last_update_time = current_time
        
        # Final summary if complete
        if self.current_step == self.total_steps:
            self._log_completion(elapsed)
    
    def _log_status(self, message: str):
        """Log a status message."""
        logging.info(message)
        if self.tty_output:
            print(f"\033[1m{message}\033[0m", file=sys.stderr)
    
    def _log_completion(self, elapsed: float):
        """Log completion message."""
        completion_msg = f"Analysis complete for {self.experiment_name} in {elapsed:.1f} seconds"
        logging.info(completion_msg)
        
        if self.tty_output:
            steps_summary = ", ".join(self.completed_steps)
            print(f"\n\033[1;32mAnalysis pipeline complete!\033[0m", file=sys.stderr)
            print(f"\033[1mTotal time: {elapsed:.1f} seconds\033[0m", file=sys.stderr)
            print(f"Completed steps: {steps_summary}", file=sys.stderr)
            print("-" * 80, file=sys.stderr)
    
    def estimate_remaining_time(self) -> float:
        """
        Estimate remaining time based on average time per step.
        
        Returns:
            Estimated remaining time in seconds
        """
        if self.current_step == 0:
            return 0.0
            
        elapsed = time.time() - self.start_time
        avg_time_per_step = elapsed / self.current_step
        remaining_steps = self.total_steps - self.current_step
        
        return avg_time_per_step * remaining_steps


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
        try:
            logging.info(f"Docker SDK version: {docker.__version__}")
        except AttributeError:
            logging.warning("Docker SDK installed but version information not available")
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