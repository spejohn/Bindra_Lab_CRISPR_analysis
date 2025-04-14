#!/usr/bin/env python3
"""
Run CRISPR analysis pipeline using Snakemake workflow manager.

This script provides a simplified interface to run the Snakemake workflow
for CRISPR screen analysis, with configurable options.
"""

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
import shlex # Import shlex for safer command splitting

def setup_logging():
    """Configure basic logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler()
        ]
    )

def get_available_cores():
    """
    Determine the number of available CPU cores.
    
    Order of precedence:
    1. CRISPR_MAX_CORES environment variable
    2. System CPU count
    3. Default to 1 if unable to determine
    
    Returns:
        int: Number of available cores
    """
    # Check for environment variable
    if "CRISPR_MAX_CORES" in os.environ:
        try:
            cores = int(os.environ["CRISPR_MAX_CORES"])
            if cores > 0:
                return cores
        except ValueError:
            logging.warning(f"Invalid CRISPR_MAX_CORES value: {os.environ['CRISPR_MAX_CORES']}, using system detection")
    
    # Check system CPU count
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass
    
    # Default fallback
    return 1

def parse_args():
    """Parse command line arguments."""
    # Get the default core count for help text
    default_cores = get_available_cores()
    
    parser = argparse.ArgumentParser(description="Run CRISPR analysis pipeline using Snakemake")
    
    # Input/output arguments
    parser.add_argument("base_dir", help="Base directory containing experiment subdirectories") # Changed from input_dir for clarity
    parser.add_argument("-o", "--output-dir", help="Directory for results (default: 'crispr_analysis_pipeline_results' in parent of base_dir)")
    
    # Workflow control
    parser.add_argument("--target_experiments", nargs='+', help="List of specific experiment directory names within base_dir to process (default: all)") # New argument
    # parser.add_argument("--configfile", help="Path to Snakemake config file") # Removed - prefer dynamic config
    parser.add_argument("-j", "--cores", type=int, 
                        help=f"Number of CPU cores for Snakemake scheduler (default: auto-detected {default_cores} cores, override with -j or CRISPR_MAX_CORES env variable)")
    parser.add_argument("--profile", help="Path to Snakemake profile directory for cluster execution") # New argument
    parser.add_argument("--dryrun", action="store_true", help="Show what would be done without executing")
    
    # Analysis options (passed via --config)
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip QC analysis")
    parser.add_argument("--skip-mle", action="store_true", help="Skip MAGeCK MLE analysis")
    
    # parser.add_argument("--snakemake-args", help="Additional arguments to pass to Snakemake", default="") # Removed - prefer explicit args
    
    return parser.parse_args()

def run_snakemake(args):
    """Run the Snakemake workflow with the given arguments."""
    # Find the Snakefile
    script_dir = Path(__file__).parent
    snakefile = script_dir / "Snakefile"
    
    if not snakefile.exists():
        logging.error(f"Snakefile not found at {snakefile}")
        return 1
    
    # Set default output directory if not specified
    output_dir = args.output_dir
    if not output_dir:
        # Place results next to the base directory
        base_dir_path = Path(args.base_dir).resolve()
        # Use the original default logic
        output_dir = base_dir_path.parent / "crispr_analysis_pipeline_results"

    # Determine number of cores to use for Snakemake scheduler
    cores = args.cores if args.cores is not None else get_available_cores()
    logging.info(f"Using {cores} CPU cores for Snakemake scheduler")
    
    # Build Snakemake command parts (as a list for better handling)
    cmd_parts = [
        "snakemake",
        "-s", str(snakefile),  # Snakefile path
        "-j", str(cores),      # Number of scheduler cores/jobs (distinct from rule threads for cluster)
        "--config", # Start config definitions
        f"base_dir={args.base_dir}", # Pass base_dir
        f"output_dir={output_dir}",
        f"skip_drugz={str(args.skip_drugz).lower()}",
        f"skip_qc={str(args.skip_qc).lower()}",
        f"skip_mle={str(args.skip_mle).lower()}",
        f"use_apptainer=true" # Align with Snakefile default
    ]
    
    # Add target experiments if specified - use target_screens key for Snakefile
    if args.target_experiments:
        # Format as a Python list string for Snakemake config
        target_list_str = "[" + ",".join(f"'{exp}'" for exp in args.target_experiments) + "]"
        cmd_parts.append(f"target_screens={target_list_str}")
    else:
        # Pass None explicitly if not provided, so Snakefile default applies
         cmd_parts.append(f"target_screens=None")


    # Add profile if specified
    if args.profile:
        cmd_parts.extend(["--profile", args.profile])
        # If using profile, --cores usually refers to scheduler cores on head node, 
        # profile handles cluster core requests via {threads}/{resources}.
        # Consider removing -j or setting it low if profile manages jobs. The profile's 'jobs:' key often controls this.
        logging.info(f"Using profile: {args.profile}. Ensure profile config handles job limits and resource allocation.")
        
    # Add dryrun if specified
    if args.dryrun:
        cmd_parts.append("--dryrun")
        
    # Join the command parts into a string
    cmd = " ".join(shlex.quote(part) for part in cmd_parts) # Use shlex.quote for safety 
    
    # Run Snakemake
    logging.info(f"Running command: {cmd}")
    try:
        # Using shell=True for simplicity here, but consider direct execution with list if needed
        # Capture output for logging
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            logging.info("Snakemake finished successfully.")
            if stdout:
                 logging.info("Snakemake stdout:\n" + stdout)
            if stderr:
                 logging.warning("Snakemake stderr:\n" + stderr) # Log stderr as warning even on success
            return 0
        else:
            logging.error(f"Snakemake failed with error code {process.returncode}")
            if stdout:
                 logging.error("Snakemake stdout:\n" + stdout)
            if stderr:
                 logging.error("Snakemake stderr:\n" + stderr)
            return process.returncode

    except Exception as e:
        logging.error(f"Error running Snakemake: {str(e)}")
        return 1

def main():
    """Main entry point."""
    setup_logging()
    args = parse_args()
    
    # Check if base directory exists
    if not os.path.isdir(args.base_dir):
        logging.error(f"Base directory {args.base_dir} does not exist or is not a directory")
        return 1
        
    # Check profile directory if specified
    if args.profile and not os.path.isdir(args.profile):
         logging.error(f"Profile directory {args.profile} does not exist or is not a directory")
         return 1

    # Run Snakemake
    return run_snakemake(args)

if __name__ == "__main__":
    sys.exit(main()) 