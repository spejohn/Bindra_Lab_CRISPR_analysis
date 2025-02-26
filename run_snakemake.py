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

def setup_logging():
    """Configure basic logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler()
        ]
    )

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run CRISPR analysis pipeline using Snakemake")
    
    # Input/output arguments
    parser.add_argument("input_dir", help="Directory containing CRISPR screen data")
    parser.add_argument("-o", "--output-dir", help="Directory for results (default: input_dir/results)")
    
    # Workflow control
    parser.add_argument("--configfile", help="Path to Snakemake config file")
    parser.add_argument("-j", "--cores", type=int, default=1, help="Number of CPU cores to use")
    parser.add_argument("--dryrun", action="store_true", help="Show what would be done without executing")
    parser.add_argument("--unlock", action="store_true", help="Unlock the working directory")
    
    # Analysis options
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip QC analysis")
    parser.add_argument("--skip-mle", action="store_true", help="Skip MAGeCK MLE analysis")
    parser.add_argument("--use-docker", action="store_true", help="Use Docker containers when available")
    
    # Snakemake parameters
    parser.add_argument("--snakemake-args", help="Additional arguments to pass to Snakemake", default="")
    
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
        output_dir = os.path.join(args.input_dir, "results")
    
    # Set default config file if not specified
    configfile = args.configfile
    if not configfile:
        default_config = script_dir / "config.yaml"
        if default_config.exists():
            configfile = str(default_config)
    
    # Build Snakemake command as a single string for direct shell execution
    cmd_parts = [
        "snakemake",
        f"-s {str(snakefile)}",  # Snakefile path
        f"-j {str(args.cores)}",  # Number of cores
    ]
    
    # Add config file if specified
    if configfile:
        cmd_parts.append(f"--configfile {configfile}")
    
    # Add config parameters
    cmd_parts.extend([
        "--config",
        f"input_dir={args.input_dir}",
        f"output_dir={output_dir}",
        f"skip_drugz={str(args.skip_drugz).lower()}",
        f"skip_qc={str(args.skip_qc).lower()}",
        f"skip_mle={str(args.skip_mle).lower()}",
        f"use_docker={str(args.use_docker).lower()}"
    ])
    
    # Add dryrun if specified
    if args.dryrun:
        cmd_parts.append("--dryrun")
    
    # Add unlock if specified
    if args.unlock:
        cmd_parts.append("--unlock")
    
    # Add additional Snakemake arguments
    if args.snakemake_args:
        cmd_parts.append(args.snakemake_args)
    
    # Join the command parts
    cmd = " ".join(cmd_parts)
    
    # Run Snakemake
    logging.info(f"Running command: {cmd}")
    try:
        # Using shell=True to avoid module import issues
        result = subprocess.run(cmd, shell=True, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        logging.error(f"Snakemake failed with error code {e.returncode}")
        return e.returncode
    except Exception as e:
        logging.error(f"Error running Snakemake: {str(e)}")
        return 1

def main():
    """Main entry point."""
    setup_logging()
    args = parse_args()
    
    # Check if input directory exists
    if not os.path.exists(args.input_dir):
        logging.error(f"Input directory {args.input_dir} does not exist")
        return 1
    
    # Run Snakemake
    return run_snakemake(args)

if __name__ == "__main__":
    sys.exit(main()) 