#!/usr/bin/env python3
"""
Run CRISPR analysis pipeline on multiple screens from a base directory.

This script automatically detects screen directories and runs the pipeline on each,
tracking which screens have already been analyzed and which need processing.
"""

import os
import argparse
import logging
import sys
from typing import List, Dict, Any
import pandas as pd
from pathlib import Path

# TODO: Refactor this script. run_pipeline no longer exists. Does this call Snakemake?
# from analysis_pipeline.pipeline import run_pipeline
from core.file_handling import identify_analyzed_experiments
from core.logging_setup import setup_logging


def detect_screen_directories(base_dir: str) -> List[str]:
    """
    Detect directories that contain CRISPR screen data.
    
    Args:
        base_dir: Base directory to search for screens
        
    Returns:
        List of screen directory paths
    """
    screen_dirs = []
    
    # Look for directories containing FASTQ files or count tables
    for root, dirs, files in os.walk(base_dir):
        has_fastq = any(f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) for f in files)
        has_count = any(f.endswith(('_count.csv', '_counts.csv', '.count', '.counts')) for f in files)
        
        if has_fastq or has_count:
            rel_path = os.path.relpath(root, base_dir)
            if rel_path == '.':
                screen_dirs.append(base_dir)
            else:
                screen_dirs.append(os.path.join(base_dir, rel_path))
    
    return screen_dirs


def find_library_file(screen_dir: str) -> str:
    """
    Find a library file in the screen directory.
    
    Args:
        screen_dir: Screen directory path
        
    Returns:
        Path to library file if found, None otherwise
    """
    # Common library file patterns
    library_patterns = [
        '*library*.csv',
        '*library*.txt',
        '*guide*.csv',
        '*guide*.txt',
        '*sgRNA*.csv',
        '*sgRNA*.txt'
    ]
    
    # Check each pattern
    for pattern in library_patterns:
        matches = list(Path(screen_dir).glob(pattern))
        if matches:
            return str(matches[0])
    
    # Look in parent directory
    parent_dir = os.path.dirname(screen_dir)
    if parent_dir != screen_dir:
        for pattern in library_patterns:
            matches = list(Path(parent_dir).glob(pattern))
            if matches:
                return str(matches[0])
    
    return None


def find_contrasts_file(screen_dir: str) -> str:
    """
    Find a contrasts file in the screen directory.
    
    Args:
        screen_dir: Screen directory path
        
    Returns:
        Path to contrasts file if found, None otherwise
    """
    # Common contrasts file patterns
    contrast_patterns = [
        '*contrast*.csv',
        '*contrast*.txt',
        '*comparison*.csv',
        '*comparison*.txt',
        '*samples*.csv',
        '*samples*.txt'
    ]
    
    # Check each pattern
    for pattern in contrast_patterns:
        matches = list(Path(screen_dir).glob(pattern))
        if matches:
            return str(matches[0])
    
    # Look in parent directory
    parent_dir = os.path.dirname(screen_dir)
    if parent_dir != screen_dir:
        for pattern in contrast_patterns:
            matches = list(Path(parent_dir).glob(pattern))
            if matches:
                return str(matches[0])
    
    return None


def run_all_screens(
    base_dir: str,
    output_base_dir: str = None,
    overwrite: bool = False,
    skip_drugz: bool = False,
    skip_qc: bool = False,
    skip_mle: bool = False,
    use_docker: bool = True  # Always use Docker, kept for backward compatibility
) -> Dict[str, Any]:
    """
    Run the pipeline on all detected screens.
    
    Args:
        base_dir: Base directory containing screen directories
        output_base_dir: Base directory for output (default: same level as base_dir, named "crispr_analysis_pipeline_results")
        overwrite: Whether to overwrite existing output files
        skip_drugz: Skip DrugZ analysis
        skip_qc: Skip quality control checks
        skip_mle: Skip MAGeCK MLE analysis
        use_docker: Parameter kept for backward compatibility (always True, Docker is required)
        
    Returns:
        Dictionary of results by screen
    """
    # Set up logging
    log_file = setup_logging(base_dir, "multi_screen_analysis")
    logging.info(f"Starting multi-screen analysis from base directory: {base_dir}")
    
    # Set default output base directory
    if output_base_dir is None:
        # Get the parent directory of the base directory for the same level
        base_parent = os.path.dirname(os.path.abspath(base_dir))
        output_base_dir = os.path.join(base_parent, "crispr_analysis_pipeline_results")
    
    # Detect screen directories
    screen_dirs = detect_screen_directories(base_dir)
    logging.info(f"Detected {len(screen_dirs)} screen directories")
    
    # Check for already analyzed experiments
    analyzed = identify_analyzed_experiments(output_base_dir)
    logging.info(f"Found {len(analyzed)} previously analyzed experiments")
    
    # Initialize results dictionary
    results = {}
    screens_df = pd.DataFrame(columns=[
        "screen_dir", "library_file", "contrasts_file", "analyzed", "success"
    ])
    
    # Process each screen
    for screen_dir in screen_dirs:
        screen_name = os.path.basename(screen_dir)
        logging.info(f"Processing screen: {screen_name}")
        
        # Check if already analyzed
        if screen_name in analyzed and analyzed[screen_name]["count_files"] and not overwrite:
            logging.info(f"Screen {screen_name} already analyzed, skipping")
            screens_df = screens_df.append({
                "screen_dir": screen_dir,
                "library_file": "N/A",
                "contrasts_file": "N/A",
                "analyzed": True,
                "success": True
            }, ignore_index=True)
            continue
        
        # Find library file
        library_file = find_library_file(screen_dir)
        if not library_file:
            logging.error(f"No library file found for screen {screen_name}, skipping")
            screens_df = screens_df.append({
                "screen_dir": screen_dir,
                "library_file": None,
                "contrasts_file": None,
                "analyzed": False,
                "success": False,
                "error": "No library file found"
            }, ignore_index=True)
            continue
        
        # Find contrasts file
        contrasts_file = find_contrasts_file(screen_dir)
        if not contrasts_file:
            logging.warning(f"No contrasts file found for screen {screen_name}, differential analysis will be skipped")
        
        # Set up output directory
        output_dir = os.path.join(output_base_dir, screen_name)
        
        # Run pipeline
        try:
            # TODO: This needs to be replaced with Snakemake call or other logic
            # screen_results = run_pipeline(
            #     input_dir=screen_dir,
            #     output_dir=output_dir,
            #     library_file=library_file,
            #     experiment_name=screen_name,
            #     contrasts_file=contrasts_file,
            #     overwrite=overwrite,
            #     skip_drugz=skip_drugz,
            #     skip_qc=skip_qc,
            #     skip_mle=skip_mle,
            #     use_docker=use_docker
            # )
            screen_results = {"status": "skipped", "error": "run_pipeline function not found"} # Placeholder
            
            results[screen_name] = screen_results
            screens_df = screens_df.append({
                "screen_dir": screen_dir,
                "library_file": library_file,
                "contrasts_file": contrasts_file,
                "analyzed": True,
                "success": "error" not in screen_results
            }, ignore_index=True)
            
            logging.info(f"Completed analysis for screen {screen_name}")
            
        except Exception as e:
            logging.error(f"Error processing screen {screen_name}: {e}")
            screens_df = screens_df.append({
                "screen_dir": screen_dir,
                "library_file": library_file,
                "contrasts_file": contrasts_file,
                "analyzed": False,
                "success": False,
                "error": str(e)
            }, ignore_index=True)
    
    # Save screen processing summary
    summary_file = os.path.join(output_base_dir, "screen_processing_summary.csv")
    screens_df.to_csv(summary_file, index=False)
    logging.info(f"Saved screen processing summary to {summary_file}")
    
    # Log completion
    logging.info(f"Multi-screen analysis complete. Processed {len(screen_dirs)} screens.")
    return results


def main():
    """
    Command-line entry point for multi-screen analysis.
    """
    parser = argparse.ArgumentParser(description="Run CRISPR analysis pipeline on multiple screens")
    
    # Required argument
    parser.add_argument("base_dir", help="Base directory containing CRISPR screen directories")
    
    # Optional arguments
    parser.add_argument("-o", "--output-dir", help="Base directory for output (defaults to base_dir/results)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--skip-drugz", action="store_true", help="Skip DrugZ analysis")
    parser.add_argument("--skip-qc", action="store_true", help="Skip quality control checks")
    parser.add_argument("--skip-mle", action="store_true", help="Skip MAGeCK MLE analysis")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Set logging level
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(log_level)
    
    # Run analysis on all screens
    run_all_screens(
        base_dir=args.base_dir,
        output_base_dir=args.output_dir,
        overwrite=args.overwrite,
        skip_drugz=args.skip_drugz,
        skip_qc=args.skip_qc,
        skip_mle=args.skip_mle,
        use_docker=True
    )
    
    print("Multi-screen analysis completed")
    return 0


if __name__ == "__main__":
    sys.exit(main()) 