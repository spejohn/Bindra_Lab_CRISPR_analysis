#!/usr/bin/env python3
"""
Utility script to convert CSV-formatted input files to the tab-delimited format
required by MAGeCK and DrugZ analysis tools.

This script will:
1. Find all design matrix and contrast table files in a directory
2. Convert them from CSV to tab-delimited format
3. Save the tab-delimited files with .txt extension
4. Clean unwanted whitespace from all values to prevent validation errors

Supports two formats for contrast tables:
- Standard format: One column each for 'contrast', 'control', and 'treatment'
- Multiple column format: Multiple columns with the same name (e.g., multiple 'control' columns)
  which will be consolidated by joining their values with commas

Usage:
    python convert_input_files.py /path/to/input_directory
    python convert_input_files.py --file /path/to/specific/file.csv
"""

import os
import sys
import glob
import argparse
import logging
from pathlib import Path
from typing import List, Optional, Dict

try:
    from analysis_pipeline.core.file_handling import convert_file_to_tab_delimited
except ImportError:
    # For standalone usage
    sys.path.append(os.path.abspath(os.path.dirname(__file__)))
    from core.file_handling import convert_file_to_tab_delimited

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger()


def find_input_files(input_dir: str, experiment_name: Optional[str] = None) -> Dict[str, List[str]]:
    """
    Find design matrix and contrast table files in the input directory structure.
    
    This function handles both new experiment-specific subdirectory structure and legacy layouts.
    
    Args:
        input_dir: Path to the base input directory
        experiment_name: Optional name of a specific experiment to focus on
        
    Returns:
        Dictionary with lists of found files
    """
    base_path = Path(input_dir)
    files = {
        'design_matrices': [],
        'contrast_tables': []
    }
    
    # Define file patterns
    design_matrix_patterns = ['*design*.csv', '*design_matrix*.csv']
    contrast_table_patterns = ['*contrast*.csv', '*comparison*.csv', '*samples*.csv']
    
    # Determine directories to search
    search_dirs = []
    
    # 1. If experiment_name is specified, look in that experiment directory first
    if experiment_name:
        exp_dir = base_path / experiment_name
        if exp_dir.exists() and exp_dir.is_dir():
            logger.info(f"Found experiment directory: {exp_dir}")
            search_dirs.append(exp_dir)
        else:
            logger.warning(f"Experiment directory '{experiment_name}' not found. Checking base directory.")
            search_dirs.append(base_path)
    else:
        # 2. If no experiment_name is specified, check if there are experiment directories
        experiment_dirs = [d for d in base_path.iterdir() if d.is_dir() and d.name not in ('fastq', 'counts')]
        
        if experiment_dirs:
            logger.info(f"Found {len(experiment_dirs)} experiment directories")
            search_dirs.extend(experiment_dirs)
        
        # 3. Always include the base directory for legacy support
        search_dirs.append(base_path)
    
    # Search each directory
    for directory in search_dirs:
        logger.debug(f"Searching in directory: {directory}")
        
        # Find design matrix files
        for pattern in design_matrix_patterns:
            matches = list(directory.glob(pattern))
            if matches:
                files['design_matrices'].extend([str(f) for f in matches])
                logger.debug(f"Found {len(matches)} design matrix files in {directory}")
                
        # Find contrast table files
        for pattern in contrast_table_patterns:
            matches = list(directory.glob(pattern))
            if matches:
                files['contrast_tables'].extend([str(f) for f in matches])
                logger.debug(f"Found {len(matches)} contrast table files in {directory}")
    
    # Deduplicate in case we found the same file via multiple paths
    files['design_matrices'] = list(set(files['design_matrices']))
    files['contrast_tables'] = list(set(files['contrast_tables']))
    
    return files


def convert_files(files: List[str], output_dir: Optional[str] = None) -> List[str]:
    """
    Convert a list of files to tab-delimited format.
    
    Args:
        files: List of file paths to convert
        output_dir: Optional output directory
        
    Returns:
        List of paths to the converted files
    """
    converted_files = []
    
    for file_path in files:
        try:
            converted_file = convert_file_to_tab_delimited(file_path, output_dir)
            converted_files.append(converted_file)
            logger.info(f"Successfully converted {file_path} to {converted_file}")
        except Exception as e:
            logger.error(f"Failed to convert {file_path}: {str(e)}")
    
    return converted_files


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Convert CSV design matrices and contrast tables to tab-delimited format"
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--dir', '-d', help="Directory containing files to convert")
    input_group.add_argument('--file', '-f', help="Specific file to convert")
    parser.add_argument('--output', '-o', help="Output directory (default: same as input)")
    parser.add_argument('--experiment', '-e', help="Specific experiment to focus on within the input directory")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    
    # Support for the first argument being a directory without a flag
    if len(sys.argv) > 1 and not sys.argv[1].startswith('-') and os.path.isdir(sys.argv[1]):
        sys.argv = [sys.argv[0], '--dir'] + sys.argv[1:]
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    if args.file:
        # Convert a single file
        if not os.path.exists(args.file):
            logger.error(f"File not found: {args.file}")
            return 1
        
        if not args.file.lower().endswith('.csv'):
            logger.error(f"File is not a CSV file: {args.file}")
            return 1
        
        try:
            converted_file = convert_file_to_tab_delimited(args.file, args.output)
            logger.info(f"Successfully converted {args.file} to {converted_file}")
        except Exception as e:
            logger.error(f"Failed to convert {args.file}: {str(e)}")
            return 1
    else:
        # Convert all files in a directory
        if not os.path.exists(args.dir):
            logger.error(f"Directory not found: {args.dir}")
            return 1
        
        files = find_input_files(args.dir, args.experiment)
        
        design_count = len(files['design_matrices'])
        contrast_count = len(files['contrast_tables'])
        
        logger.info(f"Found {design_count} design matrix files and {contrast_count} contrast table files")
        
        if design_count + contrast_count == 0:
            logger.info("No files to convert")
            return 0
        
        # Convert design matrices
        if design_count > 0:
            logger.info(f"Converting {design_count} design matrix files")
            converted_designs = convert_files(files['design_matrices'], args.output)
            logger.info(f"Successfully converted {len(converted_designs)} design matrix files")
        
        # Convert contrast tables
        if contrast_count > 0:
            logger.info(f"Converting {contrast_count} contrast table files")
            converted_contrasts = convert_files(files['contrast_tables'], args.output)
            logger.info(f"Successfully converted {len(converted_contrasts)} contrast table files")
    
    logger.info("Conversion complete")
    return 0


if __name__ == "__main__":
    sys.exit(main()) 