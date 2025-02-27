# flake8: noqa
# pylint: disable=all
# type: ignore

"""
# Snakefile for CRISPR Analysis Pipeline
# This file defines the workflow rules for processing CRISPR screen data
"""

import os
import glob
import pandas as pd
from pathlib import Path

# Configuration
configfile: "config.yaml"  # Default config file - can be overridden with --configfile

# Default configuration values
config.setdefault("input_dir", "input")
config.setdefault("output_dir", "crispr_analysis_pipeline_results")
config.setdefault("experiment_name", "experiment")
config.setdefault("norm_method", "median")
config.setdefault("fdr_threshold", 0.05)

# Helper functions
def find_library_file(wildcards):
    """Find the library file for a screen"""
    screen_dir = os.path.join(config["input_dir"], wildcards.screen)
    
    # Common library file patterns
    patterns = [
        '*library*.csv', '*library*.txt',
        '*guide*.csv', '*guide*.txt',
        '*sgRNA*.csv', '*sgRNA*.txt'
    ]
    
    # Check each pattern in screen directory
    for pattern in patterns:
        matches = glob.glob(os.path.join(screen_dir, pattern))
        if matches:
            return matches[0]
    
    # Try parent directory
    for pattern in patterns:
        matches = glob.glob(os.path.join(config["input_dir"], pattern))
        if matches:
            return matches[0]
            
    # No library file found
    raise ValueError(f"No library file found for screen {wildcards.screen}")

def find_contrasts_file(wildcards):
    """Find the contrasts file for a screen"""
    screen_dir = os.path.join(config["input_dir"], wildcards.screen)
    
    # Common contrasts file patterns
    patterns = [
        '*contrast*.csv', '*contrast*.txt',
        '*comparison*.csv', '*comparison*.txt',
        '*samples*.csv', '*samples*.txt'
    ]
    
    # Check each pattern in screen directory
    for pattern in patterns:
        matches = glob.glob(os.path.join(screen_dir, pattern))
        if matches:
            return matches[0]
    
    # Try parent directory
    for pattern in patterns:
        matches = glob.glob(os.path.join(config["input_dir"], pattern))
        if matches:
            return matches[0]
            
    return None  # It's ok if no contrasts file is found

def find_design_matrix(wildcards):
    """Find the design matrix file for a screen"""
    screen_dir = os.path.join(config["input_dir"], wildcards.screen)
    
    # Common design matrix patterns
    patterns = ['*design*.txt', '*design*.csv']
    
    # Check each pattern in screen directory
    for pattern in patterns:
        matches = glob.glob(os.path.join(screen_dir, pattern))
        if matches:
            return matches[0]
    
    # Try parent directory
    for pattern in patterns:
        matches = glob.glob(os.path.join(config["input_dir"], pattern))
        if matches:
            return matches[0]
            
    return None  # It's ok if no design matrix is found

def detect_screens():
    """Detect all screen directories within the input directory"""
    screens = []
    
    # Look for directories containing FASTQ files or count tables
    for root, dirs, files in os.walk(config["input_dir"]):
        has_fastq = any(f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) for f in files)
        has_count = any(f.endswith(('_count.csv', '_counts.csv', '.count', '.counts')) for f in files)
        
        if has_fastq or has_count:
            rel_path = os.path.relpath(root, config["input_dir"])
            if rel_path == '.':
                # If the input_dir itself is a screen
                screens.append(os.path.basename(config["input_dir"]))
            else:
                screens.append(rel_path)
    
    return screens

# Determine all screens to process
SCREENS = detect_screens()

# Default target rule
rule all:
    input:
        # Count tables for all screens
        expand(os.path.join(config["output_dir"], "{screen}", "counts", "{screen}_counts.txt"), screen=SCREENS),
        # MAGeCK RRA for all screens that have contrast files
        expand(os.path.join(config["output_dir"], "{screen}", "RRA", "{screen}.gene_summary.txt"), screen=SCREENS),
        # QC reports for all screens
        expand(os.path.join(config["output_dir"], "{screen}", "qc", "{screen}_qc_report.html"), screen=SCREENS)

# Rule to create experiment directories
rule create_experiment_dirs:
    output:
        directory(os.path.join(config["output_dir"], "{screen}")),
        directory(os.path.join(config["output_dir"], "{screen}", "counts")),
        directory(os.path.join(config["output_dir"], "{screen}", "RRA")),
        directory(os.path.join(config["output_dir"], "{screen}", "MLE")),
        directory(os.path.join(config["output_dir"], "{screen}", "drugz")),
        directory(os.path.join(config["output_dir"], "{screen}", "qc")),
        directory(os.path.join(config["output_dir"], "{screen}", "logs"))
    shell:
        "mkdir -p {output}"

# Rule to process raw FASTQ files into count tables
rule process_fastq_to_counts:
    input:
        fastq_dir = lambda wildcards: os.path.join(config["input_dir"], wildcards.screen),
        library = find_library_file,
        directories = os.path.join(config["output_dir"], "{screen}", "counts")
    output:
        counts = os.path.join(config["output_dir"], "{screen}", "counts", "{screen}_counts.txt")
    log:
        os.path.join(config["output_dir"], "{screen}", "logs", "process_fastq.log")
    shell:
        """
        python -m analysis_pipeline.analysis.sample_processing.process_all_samples \
            --fastq-dir {input.fastq_dir} \
            --library-file {input.library} \
            --output-dir {input.directories} \
            --experiment-name {wildcards.screen} > {log} 2>&1
        """

# Rule to run MAGeCK RRA analysis
rule run_mageck_rra:
    input:
        count_table = os.path.join(config["output_dir"], "{screen}", "counts", "{screen}_counts.txt"),
        contrasts = find_contrasts_file,
        directories = os.path.join(config["output_dir"], "{screen}", "RRA")
    output:
        gene_summary = os.path.join(config["output_dir"], "{screen}", "RRA", "{screen}.gene_summary.txt"),
        sgrna_summary = os.path.join(config["output_dir"], "{screen}", "RRA", "{screen}.sgrna_summary.txt")
    log:
        os.path.join(config["output_dir"], "{screen}", "logs", "mageck_rra.log")
    run:
        if input.contrasts:
            shell("""
                python -m analysis_pipeline.analysis.mageck_analysis.process_contrasts \
                    --count-table {input.count_table} \
                    --contrasts-file {input.contrasts} \
                    --output-dir {input.directories} \
                    --norm-method {config[norm_method]} \
                    --experiment-name {wildcards.screen} > {log} 2>&1
            """)
        else:
            # Create empty output files if no contrasts file exists
            shell("touch {output.gene_summary} {output.sgrna_summary}")

# Rule to run MAGeCK MLE analysis if design matrix exists
rule run_mageck_mle:
    input:
        count_table = os.path.join(config["output_dir"], "{screen}", "counts", "{screen}_counts.txt"),
        design_matrix = find_design_matrix,
        directories = os.path.join(config["output_dir"], "{screen}", "MLE")
    output:
        gene_summary = os.path.join(config["output_dir"], "{screen}", "MLE", "{screen}.gene_summary.txt")
    log:
        os.path.join(config["output_dir"], "{screen}", "logs", "mageck_mle.log")
    run:
        if input.design_matrix:
            shell("""
                python -m analysis_pipeline.analysis.mageck_analysis.process_mle_analysis \
                    --count-table {input.count_table} \
                    --design-matrix {input.design_matrix} \
                    --output-dir {input.directories} \
                    --experiment-name {wildcards.screen} \
                    --norm-method {config[norm_method]} > {log} 2>&1
            """)
        else:
            # Create empty output file if no design matrix exists
            shell("touch {output.gene_summary}")

# Rule to run DrugZ analysis
rule run_drugz:
    input:
        count_table = os.path.join(config["output_dir"], "{screen}", "counts", "{screen}_counts.txt"),
        contrasts = find_contrasts_file,
        directories = os.path.join(config["output_dir"], "{screen}", "drugz")
    output:
        results = os.path.join(config["output_dir"], "{screen}", "drugz", "{screen}_drugz_results.txt")
    log:
        os.path.join(config["output_dir"], "{screen}", "logs", "drugz.log")
    run:
        if input.contrasts:
            shell("""
                python -m analysis_pipeline.analysis.mageck_analysis.run_drugz_analysis \
                    --count-table {input.count_table} \
                    --contrasts-file {input.contrasts} \
                    --output-dir {input.directories} \
                    --prefix {wildcards.screen} > {log} 2>&1
            """)
        else:
            # Create empty output file if no contrasts file exists
            shell("touch {output.results}")

# Rule to generate QC plots and report
rule generate_qc:
    input:
        count_table = os.path.join(config["output_dir"], "{screen}", "counts", "{screen}_counts.txt"),
        rra_results = os.path.join(config["output_dir"], "{screen}", "RRA", "{screen}.gene_summary.txt"),
        mle_results = os.path.join(config["output_dir"], "{screen}", "MLE", "{screen}.gene_summary.txt"),
        drugz_results = os.path.join(config["output_dir"], "{screen}", "drugz", "{screen}_drugz_results.txt"),
        directories = os.path.join(config["output_dir"], "{screen}", "qc")
    output:
        report = os.path.join(config["output_dir"], "{screen}", "qc", "{screen}_qc_report.html")
    log:
        os.path.join(config["output_dir"], "{screen}", "logs", "qc.log")
    shell:
        """
        python -m analysis_pipeline.qc.plot_quality \
            --count-file {input.count_table} \
            --rra-file {input.rra_results} \
            --mle-file {input.mle_results} \
            --drugz-file {input.drugz_results} \
            --output-dir {input.directories} \
            --output-prefix {wildcards.screen} \
            --generate-report > {log} 2>&1
        """

# Rule to generate a summary report for all screens
rule generate_summary_report:
    input:
        qc_reports = expand(os.path.join(config["output_dir"], "{screen}", "qc", "{screen}_qc_report.html"), screen=SCREENS)
    output:
        report = os.path.join(config["output_dir"], "summary_report.html")
    log:
        os.path.join(config["output_dir"], "summary_report.log")
    shell:
        """
        python -m analysis_pipeline.reporting.generate_summary \
            --output-dir {config[output_dir]} \
            --screens {SCREENS} \
            --output-file {output.report} > {log} 2>&1
        """