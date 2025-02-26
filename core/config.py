"""
Configuration settings for CRISPR analysis pipeline.
"""

import os
from pathlib import Path

# Docker settings
DOCKER_IMAGE = "spejohn/mageck"  # Docker image for MAGeCK
DOCKER_RESTART_POLICY = {"name": "no"}  # Don't restart containers

# Default parameters
DEFAULT_NORM_METHOD = "median"  # Normalization method for MAGeCK
DEFAULT_ADJUST_METHOD = "fdr"  # Multiple testing adjustment method
DEFAULT_FDR_THRESHOLD = 0.05  # FDR threshold for significant hits

# File patterns
FASTQ_PATTERNS = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
CONTRAST_TABLE_PATTERN = "*contrast*.txt"
READ_COUNT_PATTERN = "*count*.csv"
RESULTS_MAGECK_PATTERN = "*_gMGK.csv"
RESULTS_DRUGZ_PATTERN = "*_gDZ.csv"
DESIGN_MATRIX_PATTERN = "*design*.txt"  # New pattern for design matrix files

# Required columns
LIBRARY_REQUIRED_COLUMNS = ["sgRNA", "gene"]
CONTRAST_REQUIRED_COLUMNS = ["contrast", "treatment", "control"]
SAMPLE_SHEET_COLUMNS = ["sample_name", "fastq_path"]

# Sample naming patterns
SAMPLE_NAMING_PATTERNS = {
    "Standard (_R1)": r'(.+?)_R[12]',
    "With underscore (_R1_)": r'(.+?)_R[12]_',
    "With dot (.R1)": r'(.+?)\.R[12]',
    "Standard (R1)": r'(.+?)R[12]',
    "Illumina style": r'(.+?)_S\d+_L\d+_R[12]',
}

# Essential genes for QC
DEFAULT_ESSENTIAL_GENES = {
    "Hart2015": [
        "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8", "RPL9", "RPL10", "RPL11", 
        "RPL12", "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL19", "RPL21", "RPL22", 
        "RPL23", "RPL24", "RPL26", "RPL27", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", 
        "RPL36", "RPL37", "RPL38", "RPL39", "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", 
        "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", 
        "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", 
        "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28", "RPS29", "RPS30",
        "POLR1A", "POLR1B", "POLR1C", "POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E",
        "POLR2G", "POLR2H", "POLR2I", "POLR2L", "POLR3A", "POLR3B", "POLR3C", "POLR3D",
    ]
}

# Paths
def get_drugz_path():
    """Get the path to DrugZ script, searching in common locations."""
    possible_paths = [
        Path(os.getcwd()) / "drugz" / "drugz.py",
        Path(os.getcwd()).parent / "drugz" / "drugz.py",
    ]
    
    # Also search recursively from current directory
    possible_recursive = list(Path(os.getcwd()).rglob("drugz.py"))
    possible_paths.extend(possible_recursive)
    
    for path in possible_paths:
        if path.exists():
            return path
    
    return None

# Timeouts
DOCKER_RUN_TIMEOUT = 3600  # 1 hour in seconds
DRUGZ_RUN_TIMEOUT = 900    # 15 minutes in seconds

# MLE specific configuration
DEFAULT_MLE_NORM_METHOD = "median"  # Default normalization method for MLE analysis
MLE_DESIGN_MATRIX_PATTERNS = [      # Common file patterns for design matrix files
    "*design_matrix.txt",
    "*design.txt",
    "*design_matrix.csv",
    "*design.csv"
]