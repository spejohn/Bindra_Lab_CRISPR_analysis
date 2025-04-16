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
import datetime # Import datetime for timestamps

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
    
    # Containerization flags (Mutually exclusive)
    container_group = parser.add_mutually_exclusive_group()
    container_group.add_argument("--use-apptainer", action="store_true", help="Use Apptainer for container execution")
    container_group.add_argument("--use-docker", action="store_true", help="Use Docker for container execution (default if available and neither flag set)")
    
    # Analysis options (passed via --config)
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip QC analysis")
    parser.add_argument("--skip-mle", action="store_true", help="Skip MAGeCK MLE analysis")
    
    # Add argument to capture target files/rules
    parser.add_argument('targets', nargs='*', default=[], # Capture zero or more positional args
                        help='Optional target rule names or file paths for Snakemake')
    
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
    cores = args.cores if args.cores is not None else 1 # can set to node default with get_available_cores()
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
        logging.info(f"Using profile: {args.profile}")
        # The profile now implicitly activates the executor if --executor isn't also given
        # Explicitly adding --executor is clearer
        cmd_parts.extend(["--executor", "slurm"])
        logging.info(f"Using executor: slurm (via profile or explicit flag)")

    # Add containerization flag if specified
    if args.use_apptainer:
        cmd_parts.append("--use-apptainer")
        # Optionally add Apptainer specific args if needed
        # cmd_parts.extend(["--apptainer-args", "'--bind /path/on/host:/path/in/container'"]) 
    elif args.use_docker:
        cmd_parts.append("--use-docker")
    # else: Default behavior is no container flag, snakemake runs locally

    # Add dryrun if specified
    if args.dryrun:
        cmd_parts.append("--dryrun")
        
    # Add positional targets if provided
    if args.targets:
        cmd_parts.extend(args.targets)
        logging.info(f"Adding specific targets: {args.targets}")

    # Join the command parts into a string
    # cmd = " ".join(shlex.quote(part) for part in cmd_parts) # Use shlex.quote for safety 
    # Use cmd_parts directly with Popen
    
    # Run Snakemake
    logging.info(f"Running command: {' '.join(shlex.quote(part) for part in cmd_parts)}") # Log the quoted command
    try:
        # Using shell=False and passing cmd_parts list is safer
        # Redirect stderr to stdout to capture everything in one stream
        process = subprocess.Popen(cmd_parts, 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.STDOUT, 
                                   text=True, 
                                   bufsize=1) # Line buffered

        # Read and log output line by line in real-time
        with process.stdout:
            for line in iter(process.stdout.readline, ''):
                logging.info(line.strip()) # Log each line as it arrives

        # Wait for the process to complete and get the return code
        return_code = process.wait()
        
        if return_code == 0:
            logging.info("Snakemake finished successfully.")
            return 0
        else:
            logging.error(f"Snakemake failed with error code {return_code}")
            return return_code

    except FileNotFoundError: # Handle case where snakemake command is not found
        logging.error(f"Error: 'snakemake' command not found. Is Snakemake installed and in your PATH?")
        return 1
    except Exception as e:
        logging.error(f"Error running Snakemake: {str(e)}")
        return 1

def main():
    """Main entry point."""
    setup_logging()
    args = parse_args()
    
    logging.info(f"--- Starting run_snakemake.py at {datetime.datetime.now()} ---")
    
    # Check if base directory exists
    logging.info(f"Checking base directory: {args.base_dir}")
    if not os.path.isdir(args.base_dir):
        logging.error(f"Base directory {args.base_dir} does not exist or is not a directory")
        return 1
        
    # Check profile directory if specified
    if args.profile:
        logging.info(f"Checking profile directory: {args.profile}")
        logging.info(f"Profile directory exists: {os.path.isdir(args.profile)}")
        if not os.path.isdir(args.profile):
            logging.error(f"Profile directory {args.profile} does not exist or is not a directory")
            return 1

    # Run Snakemake
    logging.info("Proceeding to execute Snakemake workflow...")
    exit_code = run_snakemake(args)
    
    logging.info(f"--- run_snakemake.py finished at {datetime.datetime.now()} with exit code {exit_code} ---")
    return exit_code

if __name__ == "__main__":
    sys.exit(main()) 