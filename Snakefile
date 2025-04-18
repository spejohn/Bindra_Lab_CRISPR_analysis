# flake8: noqa
# pylint: disable=all
# type: ignore

"""
Snakefile for CRISPR Analysis Pipeline (Refactored)
"""

# --- Rule Order Directive ---
# Resolve ambiguity between rules producing {experiment}.count.txt
ruleorder: aggregate_counts > convert_read_count_input

import os
import glob
from pathlib import Path
import sys  # Added for sys.exit
import pandas as pd  # Needed for contrast parsing
import shutil  # Added for shutil.copy
from types import SimpleNamespace  # Needed for mocking wildcards
import shlex  # Added for shlex.quote
from datetime import datetime # Import datetime
import logging  # Added for logging
import coloredlogs # Optional: for colored logs
import textwrap

# --- Logging Setup ---
# Configure basic logging
# Using coloredlogs if available, otherwise basic config
log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
# Get log level from config, default to INFO
log_level_str = config.get("log_level", "INFO").upper()
log_level = getattr(logging, log_level_str, logging.INFO) # Convert string to logging level

# TODO: Consider making log level configurable via config file or command line
try:
    coloredlogs.install(level=log_level, fmt=log_format)
    logging.info(f"Colored logging enabled at level: {log_level_str}")
except ImportError:
    logging.basicConfig(level=log_level, format=log_format)
    logging.warning(f"coloredlogs package not found, using basic logging at level: {log_level_str}")

# --- Core Function Imports ---
# Assuming core modules are importable relative to Snakefile location or via PYTHONPATH
# Adjust paths as needed based on project structure
try:
    # Use absolute import path assuming Snakefile is at root
    from core.validation import validate_experiment_structure
    from core.file_handling import convert_file_to_tab_delimited, parse_contrasts, make_count_table, convert_results_to_csv, aggregate_mageck_summaries, parse_contrast_names_from_csv
    from core.analysis_steps import (
        run_fastqc,
        run_mageck_count,
        aggregate_mageck_counts,
        run_mageck_rra,
        run_mageck_mle,
        run_drugz,
        plot_sgRNA_distribution,
        plot_gene_distribution,
        plot_gini_index,
        plot_roc_curve,
        generate_qc_report,
        # Add other imports if needed
    )

    CORE_FUNCTIONS_IMPORTED = True
    logging.info("Core functions imported successfully.")
except ImportError as e:
    logging.error(f"Error importing core functions: {e}\nEnsure the Snakefile can access the 'core' module (e.g., via PYTHONPATH or package structure).")
    CORE_FUNCTIONS_IMPORTED = False
    # Exit if core functions are missing, as the refactored rules depend on them
    sys.exit(1)

# --- Configuration ---
# Config is primarily passed via --config flag by run_snakemake.py

# Set defaults for expected config values if not provided
config.setdefault(
    "base_dir", "input"
)  # Expecting the Main_dir containing experiment subdirs
config.setdefault("output_dir", "crispr_analysis_pipeline_results")
config.setdefault(
    "target_experiments", None
)  # Optional list of specific experiment names
config.setdefault("skip_qc", False)
config.setdefault("skip_rra", False) # Add skip_rra flag
config.setdefault("skip_mle", False)
config.setdefault("skip_drugz", False)
# config.setdefault("use_apptainer", True)  # REMOVED - Containerization controlled by --use-apptainer/--use-docker

# Add defaults for other analysis options (e.g., mageck norm_method) if needed
config.setdefault("mageck_norm_method", "median")
config.setdefault("mageck_mle_options", {})
config.setdefault("drugz_options", {})

# --- Container Image Definitions ---
# These can be overridden via config file or --config flag
# *** POINT TO DOCKER HUB URIs ***
config.setdefault("fastqc_docker_uri", "docker://spejohn/crispr-analysis-fastqc:latest")
config.setdefault("mageck_docker_uri", "docker://spejohn/crispr-analysis-mageck:latest")
config.setdefault("drugz_docker_uri", "docker://spejohn/crispr-analysis-drugz:latest")

# --- Local SIF File Paths ---
# Define where the SIF files will be stored locally
SIF_DIR = Path("containers") # Store SIF files in a dedicated 'containers' directory
FASTQC_SIF = SIF_DIR / "fastqc.sif"
MAGECK_SIF = SIF_DIR / "mageck.sif"
DRUGZ_SIF = SIF_DIR / "drugz.sif"

# --- Helper Variables ---
BASE_DIR = Path(config["base_dir"]).resolve()
OUTPUT_DIR = Path(config["output_dir"]).resolve()
logging.info(f"Base Input Directory: {BASE_DIR}")
logging.info(f"Base Output Directory: {OUTPUT_DIR}")

# --- Rule to Build Apptainer SIF Files from Docker Hub ---
# This rule builds the SIF files if they don't exist or if specified differently.
# TODO: Add input dependency if needed (e.g., a version file)
rule build_sif_files:
    output:
        fastqc=FASTQC_SIF,
        mageck=MAGECK_SIF,
        drugz=DRUGZ_SIF,
    log:
        # Log file placed within the SIF directory
        str(SIF_DIR / "sif_build.log"),
    params:
        fastqc_uri=config["fastqc_docker_uri"],
        mageck_uri=config["mageck_docker_uri"],
        drugz_uri=config["drugz_docker_uri"],
        sif_dir=SIF_DIR,
    shell: r"""
        # Create the SIF directory
        mkdir -p {params.sif_dir}
        # Build each SIF file
        # Use --force to overwrite existing SIF files if the rule is triggered
        echo 'Building FastQC SIF...' >> {log}
        apptainer build --force {output.fastqc} {params.fastqc_uri} >> {log} 2>&1 && \
        echo 'Building MAGeCK SIF...' >> {log}
        apptainer build --force {output.mageck} {params.mageck_uri} >> {log} 2>&1 && \
        echo 'Building DrugZ SIF...' >> {log}
        apptainer build --force {output.drugz} {params.drugz_uri} >> {log} 2>&1
    """

# --- Experiment Discovery & Filtering ---


def find_experiments(base_dir):
    """Find potential experiment directories (immediate subdirectories)."""
    try:
        return [d.name for d in Path(base_dir).iterdir() if d.is_dir()]
    except FileNotFoundError:
        logging.error(f"Base directory not found: {base_dir}")
        sys.exit(1)


ALL_POTENTIAL_EXPERIMENTS = find_experiments(BASE_DIR)
logging.info(f"Found potential experiments: {ALL_POTENTIAL_EXPERIMENTS}")

# Filter experiments if target_screens is specified
if config["target_screens"] is not None:
    if isinstance(config["target_screens"], list):
        TARGET_EXPERIMENTS = [
            exp for exp in ALL_POTENTIAL_EXPERIMENTS if exp in config["target_screens"]
        ]
        logging.info(f"Filtering for target screens: {config['target_screens']}")
    else:
        logging.warning("target_screens in config is not a list, ignoring filter.")
        TARGET_EXPERIMENTS = ALL_POTENTIAL_EXPERIMENTS
else:
    TARGET_EXPERIMENTS = ALL_POTENTIAL_EXPERIMENTS

if not TARGET_EXPERIMENTS:
    logging.warning("No target experiments found or selected.")
    TARGET_EXPERIMENTS = []  # Ensure it's an empty list for expand

logging.info(f"Selected experiments for processing: {TARGET_EXPERIMENTS}")

# --- Validation & Input File Handling ---

# Store validation results globally to avoid re-running
ValidationCache = {}


def get_validation_info(experiment):
    global ValidationCache
    # Check for slash in experiment name
    if '/' in experiment or '\\\\' in experiment:
        # Use logging.critical for severe warnings that might cause failure
        logging.critical(f"get_validation_info called with potentially invalid experiment name containing slash: '{experiment}'")
        # Return a failure state immediately to prevent downstream issues.
        return {"status": "failed", "error": f"Invalid experiment name format: '{experiment}'"}

    if experiment not in ValidationCache:
        logging.debug(f"Cache miss for experiment: '{experiment}'") # Added quotes for clarity
        try:
            exp_path = str(BASE_DIR / experiment)
            logging.debug(f"Validating structure for path: {exp_path}") # Log path being checked
            info = validate_experiment_structure(exp_path) # Assuming this function is robust
            ValidationCache[experiment] = info
            logging.debug(f"Caching validation result for '{experiment}': {info}") # Log cached result
        except (ValueError, FileNotFoundError) as e:
            logging.error(f"Validation failed for '{experiment}': {e}")
            ValidationCache[experiment] = {"status": "failed", "error": str(e)}
            logging.debug(f"Caching validation failure for '{experiment}'") # Log failure cache
    else:
         logging.debug(f"Cache hit for experiment: '{experiment}'")

    # Log the value being returned
    result = ValidationCache.get(experiment)
    logging.debug(f"Returning validation info for '{experiment}': {result}")
    return result


# --- Helper function to dynamically determine the converted contrast/design matrix path ---
# REMOVED - No longer needed with static output paths
# def get_converted_metadata_path(wildcards, metadata_type): ... 


# --- Helper to get FASTQ basenames for an experiment ---
def get_fastq_basenames(experiment):
    """Finds sample basenames from _R1_ FASTQ files in the experiment's fastq dir."""
    validation_info = get_validation_info(experiment)
    if (
        validation_info["status"] == "failed"
        or validation_info.get("data_type") != "fastq"
    ):
        return []  # Return empty list if validation failed or not FASTQ type

    fastq_dir = validation_info.get("fastq_dir")
    if not fastq_dir or not Path(fastq_dir).is_dir():
        logging.warning(
            f"FASTQ directory not found or invalid for {experiment}: {fastq_dir}"
        )
        return []

    # Find R1 files and extract base names
    # Assumes naming convention like SampleA_R1_001.fastq.gz
    r1_files = glob.glob(str(Path(fastq_dir) / "*_R1_*.fastq.gz"))
    r1_files.extend(glob.glob(str(Path(fastq_dir) / "*_R1_*.fastq")))
    r1_files.extend(glob.glob(str(Path(fastq_dir) / "*_R1.fastq.gz")))
    r1_files.extend(glob.glob(str(Path(fastq_dir) / "*_R1.fastq")))  # Add non-001 cases

    basenames = set()
    for f in r1_files:
        # Extract the part before _R1
        basename = Path(f).name.split("_R1")[0]
        basenames.add(basename)

    logging.debug(f"Found FASTQ basenames for {experiment}: {list(basenames)}")
    return list(basenames)


# --- Rule to Convert Contrast CSV to TXT (NOW A REGULAR RULE) ---
# checkpoint convert_contrasts:
rule convert_contrasts:
    input:
        # Use the path identified by the validation function
        csv=lambda wc: get_validation_info(wc.experiment).get("contrasts_path"),
    output:
        # Static output path
        txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
    log:
        # Static log path
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_contrasts.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)
        if validation_info["status"] == "failed":
            error_msg = f"Skipping contrast conversion for {wildcards.experiment} due to validation failure: {validation_info.get('error')}"
            logging.warning(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use warning level
            with open(log[0], "w") as f:
                f.write(error_msg)
            # Exit cleanly if validation failed
            sys.exit(0)
        else:
            logging.info(f"[{datetime.now()}] {wildcards.experiment}: Converting contrasts...")
            # Output dir is the experiment dir
            output_dir_path = Path(output.txt).parent # Get directory from output file
            output_dir_path.mkdir(parents=True, exist_ok=True)
            try:
                input_csv_path = validation_info["contrasts_path"]
                if Path(input_csv_path).suffix.lower() == ".csv":
                    # Pass the FULL static output path
                    convert_file_to_tab_delimited(
                        file_path=input_csv_path,
                        output_path=str(output.txt), # Pass static path
                    )
                else:
                    logging.info( # Use info level
                        f"[{datetime.now()}] {wildcards.experiment}: Input contrast file {input_csv_path} is already .txt, copying to {output.txt}."
                    )
                    # Ensure target dir exists before copy
                    Path(output.txt).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy(input_csv_path, output.txt)
            except Exception as e:
                logging.error(f"[{datetime.now()}] converting contrasts for {wildcards.experiment}: {e}", exc_info=True) # Log exception info
                with open(log[0], "w") as f:
                    f.write(f"Error: {e}")
                raise e # Re-raise exception to fail the job


# --- Rule to Convert Design Matrix CSV to TXT (Conditional) ---
rule convert_design_matrix:
    input:
        # Use the path identified by the validation function
        csv=lambda wc: get_validation_info(wc.experiment).get("design_matrix_path")
    output:
        # Static output path
        txt=OUTPUT_DIR / "{experiment}" / "design_matrix.txt",
    # params: # Removed - no longer needed for static log path
    #     log_stem=lambda wc: Path(get_validation_info(wc.experiment).get('design_matrix_path', 'NO_INPUT')).stem
    log:
        # Static log path
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_design_matrix.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)
        if validation_info["status"] == "failed":
            # Validation failed for the experiment, just log and keep placeholder
            logging.warning( # Use warning level
                f"[{datetime.now()}] {wildcards.experiment}: Skipping design matrix conversion due to validation failure."
            )
            sys.exit(0) # Exit cleanly if validation failed
        elif validation_info.get("design_matrix_path"):
            # Design matrix exists, proceed with conversion
            input_path = validation_info["design_matrix_path"]
            output_dir_path = Path(output.txt).parent # Get directory from output file
            output_dir_path.mkdir(parents=True, exist_ok=True)
            if Path(input_path).suffix.lower() == ".csv":
                logging.info(f"[{datetime.now()}] {wildcards.experiment}: Converting design matrix...")
                try:
                    # Pass the FULL static output path
                    convert_file_to_tab_delimited(
                        file_path=input_path,
                        output_path=str(output.txt), # Pass static path
                    )
                except Exception as e:
                    logging.error( # Use error level
                        f"[{datetime.now()}] converting design matrix for {wildcards.experiment}: {e}", exc_info=True # Log exception info
                    )
                    with open(log[0], "w") as f:
                        f.write(f"Error: {e}")
                    raise e # Re-raise exception to fail the job
            else:
                logging.info( # Use info level
                    f"[{datetime.now()}] {wildcards.experiment}: Input design matrix file {input_path} is already .txt, copying to {output.txt}."
                )
                # Ensure target dir exists before copy
                Path(output.txt).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(input_path, output.txt)
        else:
            logging.info( # Use info level
                f"[{datetime.now()}] {wildcards.experiment}: Design matrix not found, skipping conversion."
            )
            sys.exit(0) # Exit cleanly if no design matrix found


# --- Helper function for conditional RC input ---
def get_rc_path_or_empty(wildcards):
    """Returns the rc_path if data_type is 'rc', otherwise an empty list."""
    validation_info = get_validation_info(wildcards.experiment)
    if validation_info.get("data_type") == "rc":
        # Return the path (which might be None if validation failed for rc type)
        return validation_info.get("rc_path")
    else:
        # Return empty list for non-rc types to avoid InputFunctionException
        return []


# --- Rule to Convert Input Read Count CSV to TSV ---
rule convert_read_count_input:
    input:
        # Use helper function to handle non-rc types gracefully during DAG build
        rc_csv=get_rc_path_or_empty,
    output:
        # Output is still defined/touched by Snakemake initially
        count_txt=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_read_count.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)

        # Run block logic remains the same - handles skipping execution
        if (
            validation_info["status"] == "valid"
            and validation_info.get("data_type") == "rc"
        ):
            input_rc_path = validation_info.get("rc_path") # Re-fetch path here for clarity
            if input_rc_path and Path(input_rc_path).exists():
                logging.info( # Use info level
                    f"[{datetime.now()}] {wildcards.experiment}: Converting input read count file: {input_rc_path} to {output.count_txt}"
                )
                Path(output.count_txt).parent.mkdir(parents=True, exist_ok=True)
                try:
                    make_count_table(
                        rc_file=str(input_rc_path),
                        output_dir=str(Path(output.count_txt).parent),
                    )
                    base_name = Path(input_rc_path).stem
                    created_file_path = (
                        Path(output.count_txt).parent / f"{base_name}.count.txt"
                    )
                    target_file_path = Path(output.count_txt)

                    if (
                        created_file_path.exists()
                        and created_file_path.resolve() != target_file_path.resolve()
                    ):
                        logging.info(f"[{datetime.now()}] {wildcards.experiment}: Renaming {created_file_path} to {target_file_path}") # Use info level
                        created_file_path.rename(target_file_path)
                    elif not target_file_path.exists():
                        if not created_file_path.exists():
                            # Let the error propagate if file not created
                            raise FileNotFoundError(
                                f"make_count_table did not produce expected file: {created_file_path} or {target_file_path}"
                            )
                    logging.info(f"[{datetime.now()}] {wildcards.experiment}: Read count conversion complete: {output.count_txt}") # Use info level
                except Exception as e:
                    # Log error and raise
                    error_msg = f"Error converting input read count for {wildcards.experiment}: {e}"
                    logging.error(f"[{datetime.now()}]: {error_msg}", exc_info=True) # Use error level, log exception
                    with open(str(log), "w") as f:
                        f.write(error_msg)
                    raise e
            else:
                # Log warning and exit cleanly (no output expected)
                # This case handles if type is 'rc' but file doesn't exist
                warn_msg = f"Warning: data_type is 'rc' but rc_path not found/valid for {wildcards.experiment}. Skipping conversion."
                logging.warning(f"[{datetime.now()}] {wildcards.experiment}: {warn_msg}") # Use warning level
                with open(str(log), "w") as f:
                    f.write(warn_msg)
                # No touch() needed here, exit cleanly
                sys.exit(0)
        else:
            # Log info and exit cleanly (no output expected) - handles non 'rc' types
            info_msg = f"Skipping read count conversion for {wildcards.experiment} (data type: {validation_info.get('data_type')}, status: {validation_info.get('status')})."
            logging.info(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}") # Use info level
            with open(str(log), "w") as f:
                f.write(info_msg)
            # No touch() needed here, exit cleanly
            sys.exit(0)


# --- Function to parse contrasts from the *converted* TXT file ---
ContrastCache = {}


# Modify signature: accept input CSV path
# def get_contrast_names(wildcards, contrast_txt_path: str):
def get_contrast_names(wildcards, contrast_csv_path: str):
    """Gets contrast names by parsing the INPUT contrast CSV file."""
    global ContrastCache
    experiment = wildcards.experiment
    # Use the INPUT path as the cache key basis
    cache_key = f"{experiment}_{contrast_csv_path}" # Use CSV path in key
    if cache_key in ContrastCache:
        # Return cached names (list of strings)
        cached_result = ContrastCache[cache_key]
        if isinstance(cached_result, list):
             # Check if it was an error marker (a list containing a single string with 'error')
            if not (len(cached_result) == 1 and isinstance(cached_result[0], str) and 'error' in cached_result[0].lower()):
                return cached_result # Return valid cached list of names
        # Fall through to re-evaluate if cache contained error or invalid type
        
    # Get validation info (needed for context, though path is passed directly)
    validation_info = get_validation_info(experiment)
    if validation_info["status"] == "failed":
        ContrastCache[cache_key] = ["error_validation_failed_upstream"] # Cache error marker
        return ["error_validation_failed_upstream"] 

    # Get the INPUT csv path from validation info to pass to the parsing function
    input_csv_path = validation_info.get("contrasts_path")
    if not input_csv_path:
        logging.error(f"Contrast CSV path not found in validation info for {experiment}") # Use error level
        ContrastCache[cache_key] = ["error_no_csv_path_in_validation"]
        return ["error_no_csv_path_in_validation"]

    try:
        # Import the new CSV parsing function
        # Ensure this import is resolvable
        from core.file_handling import parse_contrast_names_from_csv 
        
        # Call the function to get names from the INPUT CSV
        # names = parse_contrast_names_from_csv(contrast_csv_path)
        names = parse_contrast_names_from_csv(input_csv_path) # Pass the input path
        
        ContrastCache[cache_key] = names # Cache the list of names (or empty list)
        return names
    except ImportError:
        logging.error(f"Could not import parse_contrast_names_from_csv from core.file_handling") # Use error level
        ContrastCache[cache_key] = ["error_import_failed"]
        return ["error_import_failed"]
    except FileNotFoundError:
        # This case should ideally be caught by validation earlier
        logging.error(f"Input contrast CSV not found by parse_contrast_names_from_csv: {contrast_csv_path}") # Use error level
        ContrastCache[cache_key] = ["error_csv_not_found"]
        return ["error_csv_not_found"]
    except Exception as e:
        # Catch other parsing errors from the new function
        logging.error(f"Error parsing contrast names from CSV {contrast_csv_path} for {experiment}: {e}", exc_info=True) # Use error level, log exception
        ContrastCache[cache_key] = ["error_parsing_csv"] # Cache error marker
        return ["error_parsing_csv"]


# --- Function to get the primary FASTQ file for a sample basename ---
def get_fastq_for_sample(wildcards):
    """Find the R1 or single-end FASTQ for a given sample basename.
    Returns empty list if experiment type is not 'fastq' or validation failed.
    """
    validation_info = get_validation_info(wildcards.experiment)
    if (
        validation_info["status"] != "valid"
        or validation_info.get("data_type") != "fastq"
    ):
        # Return empty list instead of None
        return []

    fastq_files = validation_info.get("fastq_files", [])
    # Prefer R1
    for f in fastq_files:
        if Path(f).name.startswith(wildcards.sample) and (
            "_R1" in Path(f).name or "_R2" not in Path(f).name
        ):
            return f
    # Fallback to any matching file (single-end)
    for f in fastq_files:
        if Path(f).name.startswith(wildcards.sample):
            return f
    # Return empty list if not found (though this case is unlikely if basenames were derived correctly)
    return []


# --- Function to get the R2 FASTQ file for a sample basename ---
def get_fastq_r2_for_sample(wildcards):
    """Find the R2 FASTQ for a given sample basename, if it exists.
    Returns empty list if experiment type is not 'fastq' or R2 not found.
    """
    validation_info = get_validation_info(wildcards.experiment)
    if (
        validation_info["status"] != "valid"
        or validation_info.get("data_type") != "fastq"
    ):
        # Return empty list instead of None
        return []

    fastq_files = validation_info.get("fastq_files", [])
    # Look for R2 specifically
    for f in fastq_files:
        if Path(f).name.startswith(wildcards.sample) and "_R2" in Path(f).name:
            return f
    # Return empty list if R2 not found
    return []


# --- Rule to run FastQC on individual FASTQ files ---
rule run_fastqc_per_sample:
    input:
        fastq=get_fastq_for_sample,
        sif=FASTQC_SIF, # Depend on the specific SIF file
    output:
        # Define both outputs as FastQC generates them
        html=OUTPUT_DIR / "{experiment}" / "qc" / "{sample}_fastqc.html",
        zip=OUTPUT_DIR / "{experiment}" / "qc" / "{sample}_fastqc.zip",
    threads: 4
    resources:
        mem_mb=24000,
        time_min=120
    # No params block needed
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "fastqc_{sample}.log",
    container:
        # Directly reference the SIF path variable, converted to string
        str(FASTQC_SIF)
    # Note: Snakemake automatically translates {input.*} and {output.*} paths for the container
    #       and runs the shell command in the mounted output directory.
    shell: r"""
        mkdir -p $(dirname {log}) && \
        fastqc \
            --threads {threads} \
            -o . \
            {input.fastq} \
            > {log} 2>&1
        """


# --- Rule to run MAGeCK count per sample (from FASTQ) ---
rule run_mageck_count_per_sample:
    input:
        r1=get_fastq_for_sample,
        r2=get_fastq_r2_for_sample, # May return None if single-end
        library=lambda wc: get_validation_info(wc.experiment).get("library_path"),
        sif=MAGECK_SIF, # Keep SIF as input for dependency
    output:
        # Output files go into a dedicated subdirectory
        count_txt=OUTPUT_DIR / "{experiment}" / "mageck_count_outputs" / "{sample}.count.txt",
        # Make summary persistent for aggregation
        summary=OUTPUT_DIR / "{experiment}" / "mageck_count_outputs" / "{sample}.countsummary.txt",
    params:
        # Define paths relative to container mounts (can be helpful for debugging/reference)
        r1_container_path=lambda wc, input: f"/data/fastq/{Path(str(input.r1)).name}",
        r2_container_path=lambda wc, input: f"/data/fastq/{Path(str(input.r2)).name}" if input.r2 else "",
        library_container_path=lambda wc, input: f"/data/library/{Path(str(input.library)).name}",
        # Define host paths for binding (can be helpful for debugging/reference)
        fastq_dir_host=lambda wc, input: str(Path(str(input.r1)).parent),
        library_dir_host=lambda wc, input: str(Path(str(input.library)).parent),
        output_dir_host=lambda wc, output: str(Path(str(output.count_txt)).parent),
        # Calculate absolute output prefix string in params - Construct from output dir and sample wildcard
        output_prefix_abs=lambda wc, output: str(Path(output.count_txt).parent / wc.sample),
        # Mageck options from config
        count_options_str=lambda wc: format_options(config.get("mageck_count_options", {})),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_count_{sample}.log",
    threads: 1 # MAGeCK count is typically single-threaded
    resources:
        mem_mb=8000, # Adjust as needed
        time_min=120 # Adjust as needed
    container:
        # Directly reference the SIF path variable, converted to string
        str(MAGECK_SIF)
    # Note: Snakemake automatically translates paths and runs the shell command
    #       in the mounted host output directory (params.output_dir_host).
    #       The --output-prefix uses a relative path within that directory.
    #       Removed the conditional check for R2 input.
    shell: r"""
        mkdir -p $(dirname {log}) && \
        mageck count \
            --fastq {input.r1} \
            --list-seq {input.library} \
            --sample-label {wildcards.sample} \
            --output-prefix {params.output_prefix_abs} \
            > {log} 2>&1
        """


# --- Rule to Aggregate MAGeCK Count Summaries (within mageck_count_outputs/) ---
rule aggregate_count_summaries:
    input:
        # Input depends on the persistent summary files from run_mageck_count_per_sample
        summaries=lambda wc: expand(
            OUTPUT_DIR / "{experiment}" / "mageck_count_outputs" / "{sample}.countsummary.txt",
            experiment=wc.experiment,
            sample=get_fastq_basenames(wc.experiment),
        ) # Remove trailing comma
    output:
        # Aggregated output stays within the mageck_count_outputs directory
        agg_summary=OUTPUT_DIR / "{experiment}" / "mageck_count_outputs" / "{experiment}_aggregated.countsummary.txt",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "aggregate_count_summaries.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)

        # Only run if input type was FASTQ and validation passed
        if (
            validation_info["status"] == "valid"
            and validation_info.get("data_type") == "fastq"
        ):
            summary_files = [str(f) for f in input.summaries]
            if not summary_files:
                logging.warning( # Use warning level
                    f"[{datetime.now()}] {wildcards.experiment}: Warning: No sample count summary files found to aggregate in mageck_count_outputs/."
                )
                # Create an empty output file if needed by potential downstream rules
                Path(output.agg_summary).parent.mkdir(parents=True, exist_ok=True)
                Path(output.agg_summary).touch()
            else:
                logging.info( # Use info level
                    f"[{datetime.now()}] {wildcards.experiment}: Aggregating {len(summary_files)} sample count summaries into {output.agg_summary}..."
                )
                try:
                    # Ensure the core function exists and is imported
                    from core.file_handling import aggregate_mageck_summaries

                    Path(output.agg_summary).parent.mkdir(parents=True, exist_ok=True)
                    # Call the core function
                    success, msg = aggregate_mageck_summaries(
                        summary_files=summary_files,
                        output_path=str(output.agg_summary),
                    )
                    if not success:
                        raise RuntimeError(f"Count summary aggregation failed: {msg}")
                    
                    logging.info(f"[{datetime.now()}] {wildcards.experiment}: Count summary aggregation complete: {output.agg_summary}") # Use info level

                except ImportError:
                     error_msg = "ERROR: aggregate_mageck_summaries function not found/imported from core.file_handling."
                     logging.critical(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use critical level
                     with open(str(log), "w") as f: f.write(error_msg)
                     raise NotImplementedError(error_msg)
                except Exception as e:
                    logging.error( # Use error level
                        f"[{datetime.now()}] during count summary aggregation for {wildcards.experiment}: {e}", exc_info=True # Log exception
                    )
                    with open(str(log), "w") as f:
                        f.write(f"Error: {e}")
                    raise e
        else:
             logging.info( # Use info level
                f"[{datetime.now()}] {wildcards.experiment}: Skipping count summary aggregation (data type: {validation_info.get('data_type')}, status: {validation_info.get('status')})."
            )
             # Create empty output file if skipped? Similar consideration as aggregate_counts.
             Path(output.agg_summary).parent.mkdir(parents=True, exist_ok=True)
             Path(output.agg_summary).touch()


# --- Rule to Aggregate Sample Counts (from FASTQ processing) ---
rule aggregate_counts:
    input:
        # Input now depends on the output of the new count rule in the subdirectory
        sample_counts=lambda wc: expand(
            OUTPUT_DIR / "{experiment}" / "mageck_count_outputs" / "{sample}.count.txt", # Look in subdirectory
            experiment=wc.experiment,
            sample=get_fastq_basenames(wc.experiment)
        )
    output:
        # Aggregated output is still at the experiment level
        agg_count_txt=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "aggregate_counts.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)

        if (
            validation_info["status"] == "valid"
            and validation_info.get("data_type") == "fastq"
        ):
            if not callable(aggregate_mageck_counts):
                raise NotImplementedError(
                    "aggregate_mageck_counts function is not defined or imported."
                )

            sample_files = [str(f) for f in input.sample_counts]
            if not sample_files:
                logging.warning( # Use warning level
                    f"[{datetime.now()}] {wildcards.experiment}: Warning: No sample count files found to aggregate in mageck_count_outputs/. Creating empty output."
                )
                # Ensure output directory exists even if creating empty file
                Path(output.agg_count_txt).parent.mkdir(parents=True, exist_ok=True)
                # Create an empty file to satisfy downstream rules if needed
                Path(output.agg_count_txt).touch()

            else:
                logging.info( # Use info level
                    f"[{datetime.now()}] {wildcards.experiment}: Aggregating {len(sample_files)} sample counts from mageck_count_outputs/..."
                )
                try:
                    Path(output.agg_count_txt).parent.mkdir(parents=True, exist_ok=True)
                    # Call the core function with updated output path
                    success, msg = aggregate_mageck_counts(
                        sample_count_files=sample_files,
                        output_path=str(output.agg_count_txt),
                    )
                    if not success:
                        raise RuntimeError(f"Aggregation failed: {msg}")

                    logging.info(f"[{datetime.now()}] {wildcards.experiment}: Aggregation complete: {output.agg_count_txt}") # Use info level

                    # --- REMOVED DIRECTORY DELETION ---
                    # intermediate_counts_dir = Path(sample_files[0]).parent
                    # ... removal logic removed ...

                except Exception as e:
                    logging.error( # Use error level
                        f"[{datetime.now()}] during count aggregation for {wildcards.experiment}: {e}", exc_info=True # Log exception
                    )
                    with open(str(log), "w") as f:
                        f.write(f"Error: {e}")
                    raise e
        else:
            logging.info( # Use info level
                f"[{datetime.now()}] {wildcards.experiment}: Skipping count aggregation (data type: {validation_info.get('data_type')}, status: {validation_info.get('status')})."
            )
            # Create empty output file if skipped but needed downstream? Consider implications.
            # For now, just log the skip. Downstream rules should handle missing input if appropriate.


# --- Helper Functions for Contrast/Option Formatting ---
def get_contrast_samples(wildcards, sample_type):
    """Gets comma-separated list of treatment or control samples for a contrast."""
    experiment = wildcards.experiment
    contrast_name = wildcards.contrast
    experiment_contrasts = ContrastCache.get(experiment)
    if not isinstance(experiment_contrasts, list):
        raise ValueError(f"Contrast info not loaded for experiment '{experiment}'")
    contrast_info = next((c for c in experiment_contrasts if c.get("name") == contrast_name), None)
    if not contrast_info:
        raise ValueError(f"Could not find contrast details for '{contrast_name}' in '{experiment}'")
    samples = contrast_info.get(sample_type)
    if not samples or not isinstance(samples, list):
        raise ValueError(f"Invalid or missing '{sample_type}' samples for '{contrast_name}' in '{experiment}'")
    return ",".join(samples)

def format_options(options_dict):
    """Formats a dictionary of options into a command-line string."""
    parts = []
    for key, value in options_dict.items():
        clean_key = key.lstrip('-')
        formatted_key = f"--{clean_key}"
        if value is None or value is True:
            parts.append(formatted_key)
        elif value is False:
            pass # Skip false boolean flags
        else:
            parts.append(f"{formatted_key} {shlex.quote(str(value))}") # Quote values
    return " ".join(parts)

# NEW Helper function specifically for RRA/DrugZ rules to parse the *converted* TXT
def parse_contrast_samples_from_txt(contrast_txt_path: str, contrast_name: str, sample_type: str) -> str:
    """Parses the converted contrast.txt file to get treatment/control samples."""
    try:
        df = pd.read_csv(contrast_txt_path, sep='\t')
        contrast_row = df[df['contrast'] == contrast_name] # Use 'contrast' column
        if not contrast_row.empty:
            samples_str = contrast_row.iloc[0][sample_type]
            # Basic validation: ensure it's a non-empty string
            if isinstance(samples_str, str) and samples_str.strip():
                return samples_str
            else:
                raise ValueError(f"Invalid or empty sample list found for '{sample_type}' in contrast '{contrast_name}'")
        else:
            raise ValueError(f"Contrast '{contrast_name}' not found in {contrast_txt_path}")
    except FileNotFoundError:
        raise FileNotFoundError(f"Contrast file not found: {contrast_txt_path}")
    except KeyError:
         # Raised if 'contrast', 'treatment', or 'control' columns are missing
        raise KeyError(f"Required column missing in {contrast_txt_path}. Expecting 'contrast', 'treatment', 'control'.")
    except Exception as e:
        # Catch other potential errors during parsing
        raise RuntimeError(f"Error parsing contrasts file {contrast_txt_path} for contrast '{contrast_name}': {e}") from e

# --- Rule to run MAGeCK RRA per contrast ---
rule run_mageck_rra_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
        sif=MAGECK_SIF,
    output:
        # Simplified output path - no contrast subdirectory here
        gene_summary=OUTPUT_DIR
        / "{experiment}"
        / "analysis_results"
        / "{contrast}_RRA.gene_summary.txt",
        sgrna_summary=OUTPUT_DIR
        / "{experiment}"
        / "analysis_results"
        / "{contrast}_RRA.sgrna_summary.txt",
    params:
        analysis_options_str=lambda wc: format_options(
            {**{"norm-method": config.get("mageck_norm_method", "median")},
             **config.get("mageck_rra_options", {})}),
        # REMOVED sample parsing from params
        # treatment_samples=lambda wc, input: parse_contrast_samples_from_txt(...),
        # control_samples=lambda wc, input: parse_contrast_samples_from_txt(...),
        # Output prefix remains necessary for the command string
        output_prefix_abs=lambda wc, output: str(
            Path(output.gene_summary).parent / f"{wc.contrast}_RRA"
        ),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_rra_{contrast}.log",
    threads: 1
    container:
        str(MAGECK_SIF)
    # Changed from shell to run block
    run:
        # Imports needed within the run block scope
        import pandas as pd
        import shlex
        from snakemake.shell import shell # Explicit import
        # Needed for logging inside run block
        from datetime import datetime
        from pathlib import Path

        # Parse contrasts inside the run block from input.contrasts_txt
        treatment_samples = ""
        control_samples = ""
        try:
            # Use the helper function that reads the TXT file
            treatment_samples = parse_contrast_samples_from_txt(
                contrast_txt_path=str(input.contrasts_txt),
                contrast_name=wildcards.contrast,
                sample_type='treatment'
            )
            control_samples = parse_contrast_samples_from_txt(
                contrast_txt_path=str(input.contrasts_txt),
                contrast_name=wildcards.contrast,
                sample_type='control'
            )
        except Exception as e:
             # Log error and raise
            error_msg = f"Error parsing contrasts file {input.contrasts_txt} for contrast {wildcards.contrast} during RRA run: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f: f.write(error_msg)
            raise RuntimeError(error_msg) from e # Raise exception to stop the rule

        # Construct the command string
        command = textwrap.dedent(f"""\
        mkdir -p $(dirname {log}) && mkdir -p $(dirname {output.gene_summary}) && \
        mageck test -k {input.count_file} -t {shlex.quote(treatment_samples)} -c {shlex.quote(control_samples)} -n {params.output_prefix_abs} {params.analysis_options_str} > {log} 2>&1
        """).strip()
        # Execute the command
        shell(command)


# --- Rule to run MAGeCK MLE per experiment (Conditional) ---
rule run_mageck_mle_per_experiment:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        # Depend on the static design matrix file name
        design_matrix=OUTPUT_DIR / "{experiment}" / "design_matrix.txt",
        sif=MAGECK_SIF, # Depend on the specific SIF file
    output:
        # Experiment-level outputs in analysis_results
        gene_summary=OUTPUT_DIR / "{experiment}" / "analysis_results" / "{experiment}_MLE.gene_summary.txt",
        sgrna_summary=OUTPUT_DIR / "{experiment}" / "analysis_results" / "{experiment}_MLE.sgrna_summary.txt",
        beta_coeff=OUTPUT_DIR / "{experiment}" / "analysis_results" / "{experiment}_MLE.beta_coefficients.txt",
    params:
        # Output prefix is now relative to the container's working dir (output dir)
        # output_prefix=lambda wc, output: str(
        #     Path(output.gene_summary).parent / f"{wc.experiment}_MLE"
        # ),
        # Calculate absolute output prefix string in params
        output_prefix_abs=lambda wc, output: str(Path(output.gene_summary).with_suffix('')), 
        # Format analysis options from config
        analysis_options_str=lambda wc: format_options(config.get("mageck_mle_options", {})),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_mle_{experiment}.log", # Log per experiment
    threads: 1 # MAGeCK MLE can sometimes use more, but often limited by I/O
    container:
        # Directly reference the SIF path variable, converted to string
        str(MAGECK_SIF)
    # Snakemake ensures the output directory exists and mounts it.
    # Run mageck mle using a relative output prefix.
    shell: r"""
        mkdir -p $(dirname {log}) && \
        mageck mle \
            -k {input.count_file} \
            -d {input.design_matrix} \
            -n {params.output_prefix_abs} \
            {params.analysis_options_str} \
            > {log} 2>&1
        """


# --- Rule to run DrugZ per contrast (Conditional) ---
rule run_drugz_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
        sif=DRUGZ_SIF,
    output:
        # Simplified output path - no contrast subdirectory here
        drugz_results=OUTPUT_DIR
        / "{experiment}"
        / "analysis_results"
        / "{contrast}_DrugZ.txt",
    params:
        analysis_options_str=lambda wc: format_options(config.get("drugz_options", {})),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "drugz_{contrast}.log",
    threads: 1
    container:
        str(DRUGZ_SIF)
    run:
        # Imports needed within the run block scope
        import pandas as pd
        import shlex
        from snakemake.shell import shell # Explicit import
        # Needed for logging inside run block
        from datetime import datetime
        from pathlib import Path


        # Parse contrasts inside the run block from input.contrasts_txt
        treatment_samples = ""
        control_samples = ""
        try:
            # Simple parsing assuming tab-delimited with header: name, treatment, control
            # Input contrast file name is now dynamic via the lambda function
            df = pd.read_csv(input.contrasts_txt, sep='\t') # Needs double escape for regex/string
            contrast_row = df[df['contrast'] == wildcards.contrast] # Use 'contrast' column name
            if not contrast_row.empty:
                # Assuming treatment and control columns contain comma-separated strings already
                treatment_samples = contrast_row.iloc[0]['treatment']
                control_samples = contrast_row.iloc[0]['control']
            else:
                raise ValueError(f"Contrast '{wildcards.contrast}' not found in {input.contrasts_txt}")
        except Exception as e:
             # Log error and raise
            error_msg = f"Error parsing contrasts file {input.contrasts_txt} for contrast {wildcards.contrast}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f: f.write(error_msg)
            raise RuntimeError(error_msg) from e # Raise exception to stop the rule

        # Use the absolute path for the output file
        output_filename_abs = str(output.drugz_results)
        # Get formatted options from the dictionary
        options_str = format_options(config.get("drugz_options", {}))
        # Add unpaired flag by default, unless explicitly set to False in config
        unpaired_flag = "-unpaired" if config.get("drugz_unpaired", True) else ""

        # Construct command string explicitly without relying on textwrap or implicit newlines
        command = (
            f"mkdir -p $(dirname {log}) && "
            f"mkdir -p $(dirname {output.drugz_results}) && "
            f"python /drugz/drugz.py "
            f"-i {input.count_file} "
            f"-o \"{output_filename_abs}\" "
            f"-c {shlex.quote(control_samples)} "
            f"-x {shlex.quote(treatment_samples)} "
            f"{options_str} "
            f"{unpaired_flag} "
            f"> {log} 2>&1"
        )
        # Log the command for debugging
        logging.debug(f"Executing DrugZ command: {command}")
        # Execute the command
        shell(command)


# --- Rule to Convert MAGeCK RRA Results to CSV ---
rule convert_rra_results:
    input:
        # Updated input path to match run_mageck_rra_per_contrast output
        rra_summary=OUTPUT_DIR
        / "{experiment}"
        / "analysis_results"
        / "{contrast}_RRA.gene_summary.txt",
    output:
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}_gMGK.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_rra_{contrast}.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        # If not skipped, proceed with conversion
        try:
            convert_results_to_csv(
                result_file=str(input.rra_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="RRA"
            )
            logging.info(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: Converted RRA results.") # Use info level
        except Exception as e:
            error_msg = f"Error converting RRA results for {wildcards.contrast}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}", exc_info=True) # Use error level, log exception
            # Writing to log file is usually handled by Snakemake's log directive, but keep if custom logging needed
            # with open(log.o, "a") as f:
            #     f.write(f"[{datetime.now()}] ERROR: {error_msg}\\n")
            sys.exit(1)


# --- Rule to Convert MAGeCK MLE Results to CSV (Conditional) ---
rule convert_mle_results:
    input:
        # Update input path to match run_mageck_mle_per_experiment output directory
        mle_summary=OUTPUT_DIR / "{experiment}" / "analysis_results" / "{experiment}_MLE.gene_summary.txt",
    output:
        # Place converted file directly under experiment dir
        csv_summary=OUTPUT_DIR / "{experiment}" / "{experiment}_gMLE.csv",
    log:
        # Log is per experiment
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_mle_{experiment}.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        # If not skipped, proceed with conversion
        try:
            # Correct the input variable access here as well
            convert_results_to_csv(
                result_file=str(input.mle_summary), # Use result_file arg
                output_csv_path=str(output.csv_summary), # Use output_csv_path arg
                analysis_type="MLE" # Specify analysis type
            )
            logging.info(f"[{datetime.now()}] {wildcards.experiment}: Converted MLE results.") # Use info level
        except Exception as e:
            error_msg = f"Error converting MLE results for {wildcards.experiment}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}", exc_info=True) # Use error level, log exception
            # Writing to log file is usually handled by Snakemake's log directive, but keep if custom logging needed
            # with open(log[0], "a") as f: # Assuming log is defined; might need str(log)
            #     f.write(f"[{datetime.now()}] ERROR: {error_msg}\\n")
            sys.exit(1)


# --- Rule to Convert DrugZ Results to CSV (Conditional) ---
rule convert_drugz_results:
    input:
        # Updated input path to match run_drugz_per_contrast output
        drugz_summary=OUTPUT_DIR
        / "{experiment}"
        / "analysis_results"
        / "{contrast}_DrugZ.txt",
    output:
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}_gDZ.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_drugz_{contrast}.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        if (
            not Path(input.drugz_summary).exists()
            or Path(input.drugz_summary).stat().st_size == 0
        ):
            # Log error and raise (input required)
            error_msg = f"Input DrugZ summary {input.drugz_summary} not found or empty. Cannot convert."
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}") # Use error level
            with open(str(log), "w") as f:
                f.write(error_msg)
            # sys.exit(0) # Removed redundant exit
            raise FileNotFoundError(error_msg)

        logging.info( # Use info level
            f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: Converting DrugZ results..."
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            convert_results_to_csv(
                # Use result_file instead of input_path
                result_file=str(input.drugz_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="drugz",
            )
            logging.info(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: DrugZ results converted to CSV: {output.csv_summary}") # Use info level
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_results_to_csv function not found/imported."
            )
            logging.critical(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}") # Use critical level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting DrugZ results for {wildcards.contrast}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise e


# --- QC Rules (Conditional) ---


# Rule to plot sgRNA read count distribution
rule plot_sgrna_distribution:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
    output:
        html_plot=OUTPUT_DIR
        / "{experiment}"
        / "qc"
        / "{experiment}_sgrna_distribution.html",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "plot_sgrna_distribution.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Skipping sgRNA distribution plot."
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use error level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise FileNotFoundError(error_msg)

        logging.info(f"[{datetime.now()}] {wildcards.experiment}: Plotting sgRNA distribution...") # Use info level
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_sgRNA_distribution(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting sgRNA distribution failed: {msg}"
                logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use error level
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            logging.info(f"[{datetime.now()}] {wildcards.experiment}: sgRNA distribution plot created: {output.html_plot}") # Use info level

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_sgRNA_distribution function not found/imported."
            logging.critical(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use critical level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Error plotting sgRNA distribution for {wildcards.experiment}: {e}"
            )
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise e


# Rule to plot gene-level read count distribution
rule plot_gene_distribution:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
    output:
        html_plot=OUTPUT_DIR
        / "{experiment}"
        / "qc"
        / "{experiment}_gene_distribution.html",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "plot_gene_distribution.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Skipping gene distribution plot."
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use error level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise FileNotFoundError(error_msg)

        logging.info(f"[{datetime.now()}] {wildcards.experiment}: Plotting gene distribution...") # Use info level
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_gene_distribution(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting gene distribution failed: {msg}"
                logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use error level
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            logging.info(f"[{datetime.now()}] {wildcards.experiment}: Gene distribution plot created: {output.html_plot}") # Use info level

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_gene_distribution function not found/imported."
            logging.critical(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use critical level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Error plotting gene distribution for {wildcards.experiment}: {e}"
            )
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise e


# Rule to plot Gini index / Lorenz curve from counts
rule plot_gini_index:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
    output:
        html_plot=OUTPUT_DIR / "{experiment}" / "qc" / "{experiment}_gini_index.html",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "plot_gini_index.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Skipping Gini index plot."
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use error level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise FileNotFoundError(error_msg)

        logging.info(f"[{datetime.now()}] {wildcards.experiment}: Plotting Gini index...") # Use info level
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_gini_index(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting Gini index failed: {msg}"
                logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use error level
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            logging.info(f"[{datetime.now()}] {wildcards.experiment}: Gini index plot created: {output.html_plot}") # Use info level

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_gini_index function not found/imported."
            logging.critical(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use critical level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error plotting Gini index for {wildcards.experiment}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise e


# Placeholder rule to plot ROC curve (requires controls)
rule plot_roc_curve:
    input:
        # Input should match the output of convert_rra_results
        mageck_results=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}_gMGK.csv",
        # REMOVED: known_controls is now checked in run block
        # known_controls=lambda wc: BASE_DIR / wc.experiment / "known_controls.csv",
    output:
        html_plot=OUTPUT_DIR
        / "{experiment}"
        / "qc"
        / "{experiment}_{contrast}_roc.html",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "plot_roc_{contrast}.log",
    run:
        # Define path to optional known controls file
        known_controls_path = BASE_DIR / wildcards.experiment / "known_controls.csv"

        # Check for required MAGeCK results input
        if (
            not Path(input.mageck_results).exists()
            or Path(input.mageck_results).stat().st_size == 0
        ):
            # Log error and raise (required input)
            error_msg = f"Input MAGeCK results {input.mageck_results} not found or empty. Cannot plot ROC."
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}") # Use error level
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise FileNotFoundError(error_msg)

        # Check for optional known controls file *within* the run block
        if not known_controls_path.exists():
            # Log skip and exit cleanly (optional input)
            warn_msg = f"Known controls file {known_controls_path} not found. Skipping ROC plot for {wildcards.experiment}/{wildcards.contrast}."
            logging.warning(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {warn_msg}") # Use warning level
            with open(str(log), "w") as f:
                f.write(warn_msg)
            # REMOVED: Create a dummy output file to satisfy Snakemake if skipped
            # Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            # Path(output.html_plot).touch()
            sys.exit(0) # Exit successfully without creating output

        # Proceed with plotting if controls file exists
        logging.info( # Use info level
            f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: Plotting ROC curve (PLACEHOLDER)..."
        )
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            logging.debug( # Use debug level for placeholder message
                "--> Placeholder rule: plot_roc_curve() function would be called here."
            )
            # Placeholder: Just log completion, no touch needed (dummy file created above if skipped)
            # If plotting were implemented, it would overwrite the dummy file here.
            # For now, just log completion, actual plot generation is placeholder.
            logging.info(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: ROC curve plot logged completion (placeholder): {output.html_plot}") # Use info level
            # Ensure output file exists even for placeholder logic
            Path(output.html_plot).touch()


        except NameError:
            error_msg = "ERROR: plot_roc_curve function not found/imported (placeholder assumes it exists)."
            logging.critical(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}") # Use critical level
            with open(str(log), "w") as f:
                f.write(error_msg)
            # Raise error if function is expected but missing
            raise
        except Exception as e:
            error_msg = f"Error plotting ROC curve for {wildcards.experiment}/{wildcards.contrast}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise e


# Placeholder rule to generate an aggregate QC report
rule generate_qc_report:
    input:
        sgrna_dist=OUTPUT_DIR
        / "{experiment}"
        / "qc"
        / "{experiment}_sgrna_distribution.html",
        gene_dist=OUTPUT_DIR
        / "{experiment}"
        / "qc"
        / "{experiment}_gene_distribution.html",
        gini=OUTPUT_DIR / "{experiment}" / "qc" / "{experiment}_gini_index.html",
        # Modify the lambda to fetch and pass the contrast_csv_path
        roc_curves=lambda wc: (
            # Fetch validation info and contrast path first
            validation_info := get_validation_info(wc.experiment),
            contrast_csv_path := validation_info.get("contrasts_path"),
            # Get contrast names, handle None path
            contrasts := (
                get_contrast_names(
                    SimpleNamespace(experiment=wc.experiment),
                    contrast_csv_path=contrast_csv_path
                )
                if contrast_csv_path else []
            ),
            # Perform the expand using the fetched contrasts
            expand(
                OUTPUT_DIR
                / wc.experiment
                / "qc"
                / f"{wc.experiment}_{{contrast}}_roc.html",
                contrast=contrasts,
            )
        )[-1], # Return the last element (the result of expand)
        fastqc_reports=lambda wc: (
            expand(
                OUTPUT_DIR / wc.experiment / "qc" / "{sample}_fastqc.html",
                sample=get_fastq_basenames(wc.experiment),
            )
            if get_validation_info(wc.experiment).get("data_type") == "fastq"
            else []
        ),
    output:
        report=OUTPUT_DIR / "{experiment}" / "qc" / "{experiment}_qc_report.html",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "generate_qc_report.log",
    run:
        # Removed skip check - Handled by rule all's conditional input
        logging.info(f"[{datetime.now()}] {wildcards.experiment}: Generating QC report (PLACEHOLDER)...") # Use info level
        try:
            Path(output.report).parent.mkdir(parents=True, exist_ok=True)
            input_files_dict = {
                "sgrna_dist": (
                            str(input.sgrna_dist) if Path(input.sgrna_dist).exists() else None
                        ),
                        "gene_dist": (
                            str(input.gene_dist) if Path(input.gene_dist).exists() else None
                        ),
                        "gini": str(input.gini) if Path(input.gini).exists() else None,
                        "roc_curves": [str(f) for f in input.roc_curves if Path(f).exists()],
                        "fastqc_reports": [
                            str(f) for f in input.fastqc_reports if Path(f).exists()
                        ],
            }

            logging.debug( # Use debug level for placeholder message
                "--> Placeholder rule: generate_qc_report() function would be called here with:"
            )
            logging.debug(input_files_dict) # Use debug level
            # Placeholder: Log completion, no touch needed
            logging.info(f"[{datetime.now()}] {wildcards.experiment}: QC report logged completion (placeholder): {output.report}") # Use info level

        except NameError:
            error_msg = "ERROR: generate_qc_report function not found/imported (placeholder assumes it exists)."
            logging.critical(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}") # Use critical level
            with open(str(log), "w") as f:
                f.write(error_msg)
        except Exception as e:
            error_msg = f"Error generating QC report for {wildcards.experiment}: {e}"
            logging.error(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}", exc_info=True) # Use error level, log exception
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise e


# --- Rule All: Define Final Output Files ---


def get_final_outputs(wildcards):
    """Dynamically collects all expected final output files based on config and targets."""
    final_files = []
    logging.debug("--- Entering get_final_outputs ---") # DEBUG -> debug
    # REMOVED DEBUG PRINT
    # logging.debug(f"Inside get_final_outputs: config['skip_drugz'] = {config.get('skip_drugz')}")
    logging.debug(f"Target experiments: {TARGET_EXPERIMENTS}") # DEBUG -> debug

    # Use wildcards.experiment if rule all defines experiment, otherwise iterate
    # Assuming rule all doesn't have {experiment} wildcard for simplicity now
    # We iterate over TARGET_EXPERIMENTS determined earlier
    for experiment in TARGET_EXPERIMENTS:
        logging.debug(f"Processing experiment: {experiment}") # DEBUG -> debug
        validation_info = get_validation_info(experiment)
        logging.debug(f"Validation info for {experiment}: {validation_info}") # DEBUG -> debug

        if validation_info["status"] == "failed":
            logging.warning( # Use warning level
                f"Skipping final output collection for failed experiment: {experiment}"
            )
            continue

        # Determine the path to the INPUT contrast CSV file from validation info
        # contrast_txt_path = OUTPUT_DIR / experiment / "contrasts.txt"
        input_contrast_csv_path = validation_info.get("contrasts_path")

        # Check if the input path was found during validation
        if not input_contrast_csv_path:
             logging.error(f"Could not determine input contrast CSV path for {experiment} from validation. Skipping contrast-dependent outputs.") # Use error level
             contrasts = ["error_no_csv_path_found"]
        else:
            # Use the modified get_contrast_names function
            # Pass the INPUT CSV path
            logging.debug(f"Calling get_contrast_names for {experiment} with INPUT path={input_contrast_csv_path}...") # DEBUG -> debug
            # contrasts = get_contrast_names(SimpleNamespace(experiment=experiment), contrast_txt_path=str(contrast_txt_path))
            contrasts = get_contrast_names(SimpleNamespace(experiment=experiment), contrast_csv_path=input_contrast_csv_path) # Pass CSV path
            logging.debug(f"Contrasts for {experiment}: {contrasts}") # DEBUG -> debug

        # The rest of the function uses the `contrasts` list as before
        # Error handling for contrasts list is already present

        # --- Explicitly depend on the converted contrasts.txt for this experiment ---
        # This ensures the convert_contrasts checkpoint runs before rule all input is finalized
        contrast_txt_path = OUTPUT_DIR / experiment / "contrasts.txt"
        logging.debug(f"Adding explicit dependency: {contrast_txt_path}") # DEBUG -> debug
        final_files.append(contrast_txt_path)
        # --------------------------------------------------------------------------

        if (
            not isinstance(contrasts, list) or (contrasts and "error" in contrasts[0]) # Handle error case and empty list
        ):
            logging.warning( # Use warning level
                f"Skipping contrast-specific outputs for {experiment} due to contrast parsing issues."
            )
            contrasts = [] # Ensure contrasts is an empty list if issues occurred

        # --- 1. Analysis Results (per contrast / per experiment) ---
        logging.debug(f"Checking analysis results for {experiment}...") # DEBUG -> debug
        for contrast in contrasts:
            logging.debug(f"  Checking contrast: {contrast}") # DEBUG -> debug
            # Add RRA results if not skipped
            if not config.get("skip_rra", False):
                # DEPEND ON CONVERTED CSV
                rra_file = OUTPUT_DIR / experiment / f"{contrast}_gMGK.csv"
                logging.debug(f"    Adding RRA target: {rra_file}") # DEBUG -> debug
                final_files.append(rra_file)
            else:
                logging.debug("    Skipping RRA (skip_rra=True)") # DEBUG -> debug

            # Add DrugZ results if not skipped
            # Revert check to simple boolean
            if not config.get("skip_drugz", False):
                # DEPEND ON CONVERTED CSV
                dz_file = OUTPUT_DIR / experiment / f"{contrast}_gDZ.csv"
                logging.debug(f"    Adding DrugZ target: {dz_file}") # DEBUG -> debug
                final_files.append(dz_file)
            else:
                logging.debug(f"    Skipping DrugZ for {contrast} (skip_drugz=True)") # DEBUG -> debug

        # Add MLE results (per experiment) if not skipped AND design matrix exists
        logging.debug(f"Checking MLE results for {experiment}...") # DEBUG -> debug
        # Check for *converted* design matrix in output dir
        design_matrix_path = OUTPUT_DIR / experiment / "design_matrix.txt"
        # Check if the *rule output* design matrix exists
        design_output_exists = design_matrix_path.exists()
        logging.debug(f"  MLE Check: skip_mle={config.get('skip_mle', False)}, design_matrix_exists={design_output_exists} ({design_matrix_path})") # DEBUG -> debug

        if not config.get("skip_mle", False) and design_output_exists:
             # REMOVE NATIVE MLE OUTPUTS AS TARGETS
             # mle_gene_native = OUTPUT_DIR / experiment / "analysis_results" / f"{experiment}_MLE.gene_summary.txt"
             # mle_sgrna_native = OUTPUT_DIR / experiment / "analysis_results" / f"{experiment}_MLE.sgrna_summary.txt"
             # mle_beta_native = OUTPUT_DIR / experiment / "analysis_results" / f"{experiment}_MLE.beta_coefficients.txt"
             # DEPEND ON CONVERTED MLE CSV (at experiment level)
             mle_csv = OUTPUT_DIR / experiment / f"{experiment}_gMLE.csv"
             # logging.debug(f"    Adding MLE targets: {mle_gene_native}, {mle_sgrna_native}, {mle_beta_native}, {mle_csv}") # DEBUG -> debug
             logging.debug(f"    Adding MLE target: {mle_csv}") # DEBUG -> debug
             # final_files.append(mle_gene_native)
             # final_files.append(mle_sgrna_native)
             # final_files.append(mle_beta_native)
             final_files.append(mle_csv)
        else:
            logging.debug(f"    Skipping MLE targets.") # DEBUG -> debug

        # --- 2. QC Files (Conditional) ---
        logging.debug(f"Checking QC files for {experiment}...") # DEBUG -> debug
        # Revert check to simple boolean
        if not config.get("skip_qc", False):
            # Per-experiment QC plots (ensure rules still exist)
            if "plot_sgrna_distribution" in globals(): # Check if rule exists
                sgrna_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_sgrna_distribution.html"
                logging.debug(f"    Adding QC target: {sgrna_plot}") # DEBUG -> debug
                final_files.append(sgrna_plot)
            if "plot_gene_distribution" in globals():
                gene_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_gene_distribution.html"
                logging.debug(f"    Adding QC target: {gene_plot}") # DEBUG -> debug
                final_files.append(gene_plot)
            if "plot_gini_index" in globals():
                gini_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_gini_index.html"
                logging.debug(f"    Adding QC target: {gini_plot}") # DEBUG -> debug
                final_files.append(gini_plot)
            if "generate_qc_report" in globals():
                qc_report = OUTPUT_DIR / experiment / "qc" / f"{experiment}_qc_report.html"
                logging.debug(f"    Adding QC target: {qc_report}") # DEBUG -> debug
                final_files.append(qc_report)

            # Per-contrast QC plots (ROC needs controls check)
            if "plot_roc_curve" in globals():
                 known_controls_path = BASE_DIR / experiment / "known_controls.csv"
                 logging.debug(f"  ROC Check: known_controls_exists={known_controls_path.exists()} ({known_controls_path})") # DEBUG -> debug
                 if known_controls_path.exists(): # Only *require* ROC if controls file exists
                    for contrast in contrasts:
                        # This plot is now generated conditionally by the rule itself,
                        # but we still need it as a target for rule all if controls exist.
                        roc_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_{contrast}_roc.html"
                        logging.debug(f"    Adding QC target (since controls exist): {roc_plot}") # DEBUG -> debug
                        final_files.append(roc_plot)
                 else:
                     logging.debug("    Skipping adding ROC plots to final targets (no known_controls.csv)") # DEBUG -> debug
            else:
                 logging.debug("    Skipping ROC plots (rule not defined)") # DEBUG -> debug

            # Per-sample FastQC reports (if FASTQ input)
            logging.debug(f"  FastQC Check: data_type={validation_info.get('data_type')}") # DEBUG -> debug
            if validation_info.get("data_type") == "fastq" and "run_fastqc_per_sample" in globals():
                fastq_basenames = get_fastq_basenames(experiment)
                logging.debug(f"    FastQC basenames: {fastq_basenames}") # DEBUG -> debug
                for sample in fastq_basenames:
                    fq_html = OUTPUT_DIR / experiment / "qc" / f"{sample}_fastqc.html"
                    fq_zip = OUTPUT_DIR / experiment / "qc" / f"{sample}_fastqc.zip"
                    logging.debug(f"    Adding QC targets: {fq_html}, {fq_zip}") # DEBUG -> debug
                    final_files.append(fq_html)
                    final_files.append(fq_zip)
            else:
                logging.debug("    Skipping FastQC reports (not FASTQ data or rule undefined)") # DEBUG -> debug
        else:
            logging.debug("  Skipping ALL QC (skip_qc=True)") # DEBUG -> debug

    # Convert Path objects to strings for Snakemake input list
    final_files_str = [str(f) for f in final_files]
    logging.debug(f"--- Exiting get_final_outputs with {len(final_files_str)} targets: {final_files_str} ---") # DEBUG -> debug
    return final_files_str


# --- Final Target Rule ---
# Remove the temporary get_first_experiment_qc_targets function
# def get_first_experiment_qc_targets(): ...


rule all:
    input:
        # Explicitly list outputs that depend on the converted contrasts rule
        contrast_files=expand(
            OUTPUT_DIR / "{experiment}" / "contrasts.txt",
            experiment=TARGET_EXPERIMENTS # Use the globally determined list
        ),
        # Use lambda for the rest of the outputs, remove checkpoints argument
        # other_files=lambda wildcards, checkpoints: get_final_outputs(wildcards, checkpoints),
        other_files=lambda wildcards: get_final_outputs(wildcards),
    output:
        # Single flag file indicating completion of requested targets
        OUTPUT_DIR / "pipeline_complete.flag",
    run:
        logging.info(f"[{datetime.now()}] CRISPR Analysis Pipeline Workflow Completed for targets:") # Use info level
        # Input is now structured, access specific parts if needed
        # logging.debug(f"Required contrast files: {input.contrast_files}") # Use debug level
        # logging.debug(f"Other required files: {input.other_files}") # Use debug level
        logging.info(f"Input files generated (list might be long). Check output directory.") # Use info level
        logging.info(f"Completion flag: {output}") # Use info level
        # Create the flag file
        Path(output[0]).touch()
