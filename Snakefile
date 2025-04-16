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

# --- Core Function Imports ---
# Assuming core modules are importable relative to Snakefile location or via PYTHONPATH
# Adjust paths as needed based on project structure
try:
    # Use absolute import path assuming Snakefile is at root
    from core.validation import validate_experiment_structure
    from core.file_handling import convert_file_to_tab_delimited, parse_contrasts, make_count_table, convert_results_to_csv, aggregate_mageck_summaries
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
    print("Core functions imported successfully.")
except ImportError as e:
    print(f"Error importing core functions: {e}\nEnsure the Snakefile can access the 'core' module (e.g., via PYTHONPATH or package structure).")
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
print(f"Base Input Directory: {BASE_DIR}")
print(f"Base Output Directory: {OUTPUT_DIR}")

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
        print(f"Error: Base directory not found: {base_dir}")
        sys.exit(1)


ALL_POTENTIAL_EXPERIMENTS = find_experiments(BASE_DIR)
print(f"Found potential experiments: {ALL_POTENTIAL_EXPERIMENTS}")

# Filter experiments if target_screens is specified
if config["target_screens"] is not None:
    if isinstance(config["target_screens"], list):
        TARGET_EXPERIMENTS = [
            exp for exp in ALL_POTENTIAL_EXPERIMENTS if exp in config["target_screens"]
        ]
        print(f"Filtering for target screens: {config['target_screens']}")
    else:
        print("Warning: target_screens in config is not a list, ignoring filter.")
        TARGET_EXPERIMENTS = ALL_POTENTIAL_EXPERIMENTS
else:
    TARGET_EXPERIMENTS = ALL_POTENTIAL_EXPERIMENTS

if not TARGET_EXPERIMENTS:
    print("No target experiments found or selected.")
    TARGET_EXPERIMENTS = []  # Ensure it's an empty list for expand

print(f"Selected experiments for processing: {TARGET_EXPERIMENTS}")

# --- Validation & Input File Handling ---

# Store validation results globally to avoid re-running
ValidationCache = {}


def get_validation_info(experiment):
    global ValidationCache
    if experiment not in ValidationCache:
        try:
            exp_path = str(BASE_DIR / experiment)
            info = validate_experiment_structure(exp_path)
            ValidationCache[experiment] = info
            print(f"Validation successful for {experiment}")
        except (ValueError, FileNotFoundError) as e:
            print(f"ERROR: Validation failed for {experiment}: {e}")
            ValidationCache[experiment] = {"status": "failed", "error": str(e)}
    return ValidationCache[experiment]


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
        print(
            f"Warning: FASTQ directory not found or invalid for {experiment}: {fastq_dir}"
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

    print(f"Found FASTQ basenames for {experiment}: {list(basenames)}")
    return list(basenames)


# --- Rule to Convert Contrast CSV to TXT ---
checkpoint convert_contrasts:
    input:
        # Use the path identified by the validation function
        csv=lambda wc: get_validation_info(wc.experiment).get("contrasts_path"),
    output:
        # Static output path
        txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
    # params: # Removed - no longer needed for static log path
    #     log_stem=lambda wc: Path(get_validation_info(wc.experiment).get('contrasts_path', '')).stem
    log:
        # Static log path
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_contrasts.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)
        if validation_info["status"] == "failed":
            error_msg = f"Skipping contrast conversion for {wildcards.experiment} due to validation failure: {validation_info.get('error')}"
            print(f"[{datetime.now()}] {wildcards.experiment}: {error_msg}")
            with open(log[0], "w") as f:
                f.write(error_msg)
            # Exit cleanly if validation failed
            sys.exit(0)
        else:
            print(f"[{datetime.now()}] {wildcards.experiment}: Converting contrasts...")
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
                    print(
                        f"[{datetime.now()}] {wildcards.experiment}: Input contrast file {input_csv_path} is already .txt, copying to {output.txt}."
                    )
                    # Ensure target dir exists before copy
                    Path(output.txt).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy(input_csv_path, output.txt)
            except Exception as e:
                print(f"[{datetime.now()}] ERROR converting contrasts for {wildcards.experiment}: {e}")
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
            print(
                f"[{datetime.now()}] {wildcards.experiment}: Skipping design matrix conversion due to validation failure."
            )
            sys.exit(0) # Exit cleanly if validation failed
        elif validation_info.get("design_matrix_path"):
            # Design matrix exists, proceed with conversion
            input_path = validation_info["design_matrix_path"]
            output_dir_path = Path(output.txt).parent # Get directory from output file
            output_dir_path.mkdir(parents=True, exist_ok=True)
            if Path(input_path).suffix.lower() == ".csv":
                print(f"[{datetime.now()}] {wildcards.experiment}: Converting design matrix...")
                try:
                    # Pass the FULL static output path
                    convert_file_to_tab_delimited(
                        file_path=input_path,
                        output_path=str(output.txt), # Pass static path
                    )
                except Exception as e:
                    print(
                        f"[{datetime.now()}] ERROR converting design matrix for {wildcards.experiment}: {e}"
                    )
                    with open(log[0], "w") as f:
                        f.write(f"Error: {e}")
                    raise e # Re-raise exception to fail the job
            else:
                print(
                    f"[{datetime.now()}] {wildcards.experiment}: Input design matrix file {input_path} is already .txt, copying to {output.txt}."
                )
                # Ensure target dir exists before copy
                Path(output.txt).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(input_path, output.txt)
        else:
            print(
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
                print(
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
                        print(f"[{datetime.now()}] {wildcards.experiment}: Renaming {created_file_path} to {target_file_path}")
                        created_file_path.rename(target_file_path)
                    elif not target_file_path.exists():
                        if not created_file_path.exists():
                            # Let the error propagate if file not created
                            raise FileNotFoundError(
                                f"make_count_table did not produce expected file: {created_file_path} or {target_file_path}"
                            )
                    print(f"[{datetime.now()}] {wildcards.experiment}: Read count conversion complete: {output.count_txt}")
                except Exception as e:
                    # Log error and raise
                    error_msg = f"Error converting input read count for {wildcards.experiment}: {e}"
                    print(f"[{datetime.now()}] ERROR: {error_msg}")
                    with open(str(log), "w") as f:
                        f.write(error_msg)
                    raise e
            else:
                # Log warning and exit cleanly (no output expected)
                # This case handles if type is 'rc' but file doesn't exist
                warn_msg = f"Warning: data_type is 'rc' but rc_path not found/valid for {wildcards.experiment}. Skipping conversion."
                print(f"[{datetime.now()}] {wildcards.experiment}: {warn_msg}")
                with open(str(log), "w") as f:
                    f.write(warn_msg)
                # No touch() needed here, exit cleanly
                sys.exit(0)
        else:
            # Log info and exit cleanly (no output expected) - handles non 'rc' types
            info_msg = f"Skipping read count conversion for {wildcards.experiment} (data type: {validation_info.get('data_type')}, status: {validation_info.get('status')})."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # No touch() needed here, exit cleanly
            sys.exit(0)


# --- Function to parse contrasts from the *converted* TXT file ---
ContrastCache = {}


def get_contrast_names(wildcards, checkpoints):
    """Gets contrast names by reading the output of the convert_contrasts checkpoint."""
    global ContrastCache
    experiment = wildcards.experiment
    cache_key = experiment
    if cache_key in ContrastCache:
        # Check if cached result was an error placeholder
        if isinstance(ContrastCache[cache_key], list) and ContrastCache[cache_key]:
            # Return only valid names
            return [c["name"] for c in ContrastCache[cache_key] if isinstance(c, dict)]
        elif isinstance(
            ContrastCache[cache_key], list
        ):  # Empty list means no valid contrasts
            return []
        else:  # It was an error string/placeholder
            # Rerun logic if it was a placeholder error
            pass # Fall through to re-evaluate

    # --- Access checkpoint output --- 
    try:
        # Use the checkpoints object to get the output path *after* the checkpoint runs
        converted_contrast_path = checkpoints.convert_contrasts.get(experiment=experiment).output.txt
    except Exception as e:
        # Handle cases where the checkpoint might not have been run or failed
        print(f"Error accessing checkpoint output for {experiment}: {e}")
        ContrastCache[cache_key] = ["error_accessing_checkpoint"]
        return ["error_accessing_checkpoint"]
    # --- End Checkpoint Access ---

    validation_info = get_validation_info(experiment)
    if validation_info["status"] == "failed":
        ContrastCache[cache_key] = ["validation_failed"]
        return ["validation_failed"]

    # File existence check might still be useful, but rely on checkpoint dependency primarily
    # if not Path(converted_contrast_path).exists():
    #     print(f"Warning: Converted contrasts file not found via checkpoint: {converted_contrast_path}")
    #     ContrastCache[cache_key] = ["conversion_failed"]
    #     return ["conversion_failed"]

    try:
        # Ensure the path from the checkpoint output is used
        contrasts = parse_contrasts(str(converted_contrast_path))
        ContrastCache[cache_key] = contrasts  # Cache the list of dicts
        names = [c["name"] for c in contrasts] if contrasts else []
        if not names:
            print(f"Warning: No valid contrasts parsed for {experiment}.")
            return []  # Return empty list
        return names
    except Exception as e:
        print(f"Error parsing contrasts for {experiment}: {e}")
        ContrastCache[cache_key] = ["error_parsing_contrasts"]
        return ["error_parsing_contrasts"]


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
    params: # Add params directive back, even if empty
        # No specific parameters needed for the shell command anymore
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "fastqc_{sample}.log",
    container:
        # Directly reference the SIF path variable, converted to string
        str(FASTQC_SIF)
    shell:
        # Snakemake ensures the output directory exists on the host and mounts it.
        # Run fastqc, outputting to the current directory (.) inside the container.
        r"""
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
        # Define paths relative to container mounts
        # These help define the mount points for Snakemake
        r1_container_path=lambda wc, input: f"/data/fastq/{Path(str(input.r1)).name}",
        r2_container_path=lambda wc, input: f"/data/fastq/{Path(str(input.r2)).name}" if input.r2 else "",
        library_container_path=lambda wc, input: f"/data/library/{Path(str(input.library)).name}",
        # Define the prefix path *inside* the container's output mount
        # Snakemake maps the host output dir (mageck_count_outputs) to /data/output
        output_prefix_container="/data/output/{sample}",
        # Define host paths for binding
        fastq_dir_host=lambda wc, input: str(Path(str(input.r1)).parent),
        library_dir_host=lambda wc, input: str(Path(str(input.library)).parent),
        # Host output directory is now the new subdirectory
        output_dir_host=lambda wc, output: str(Path(str(output.count_txt)).parent),
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
    shell:
        # Ensure the target output directory exists *inside* the container
        # Note: Snakemake automatically translates {input.*} and {output.*} paths for the container
        r"""
        mkdir -p $(dirname {params.output_prefix_container});
        mageck count \\
            --fastq {input.r1} \\
            $(test -n "{input.r2}" && echo "--fastq-2 {input.r2}") \\
            --list-seq {input.library} \\
            --sample-label {wildcards.sample} \\
            --output-prefix {params.output_prefix_container} \\
            {params.count_options_str} \\
            > {log} 2>&1
        # No mv commands needed now, files are created in the mounted output directory
        """


# --- Rule to Aggregate MAGeCK Count Summaries (within mageck_count_outputs/) ---
rule aggregate_count_summaries:
    input:
        # Input depends on the persistent summary files from run_mageck_count_per_sample
        summaries=lambda wc: expand(
            OUTPUT_DIR / "{experiment}" / "mageck_count_outputs" / "{sample}.countsummary.txt",
            experiment=wc.experiment,
            sample=get_fastq_basenames(wc.experiment),
        ),
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
                print(
                    f"[{datetime.now()}] {wildcards.experiment}: Warning: No sample count summary files found to aggregate in mageck_count_outputs/."                    
                )
                # Create an empty output file if needed by potential downstream rules
                Path(output.agg_summary).parent.mkdir(parents=True, exist_ok=True)
                Path(output.agg_summary).touch()
            else:
                print(
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
                    
                    print(f"[{datetime.now()}] {wildcards.experiment}: Count summary aggregation complete: {output.agg_summary}")

                except ImportError:
                     error_msg = "ERROR: aggregate_mageck_summaries function not found/imported from core.file_handling."
                     print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
                     with open(str(log), "w") as f: f.write(error_msg)
                     raise NotImplementedError(error_msg)
                except Exception as e:
                    print(
                        f"[{datetime.now()}] ERROR during count summary aggregation for {wildcards.experiment}: {e}"
                    )
                    with open(str(log), "w") as f:
                        f.write(f"Error: {e}")
                    raise e
        else:
             print(
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
            sample=get_fastq_basenames(wc.experiment),
        ),
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
                print(
                    f"[{datetime.now()}] {wildcards.experiment}: Warning: No sample count files found to aggregate in mageck_count_outputs/. Creating empty output."
                )
                # Ensure output directory exists even if creating empty file
                Path(output.agg_count_txt).parent.mkdir(parents=True, exist_ok=True)
                # Create an empty file to satisfy downstream rules if needed
                Path(output.agg_count_txt).touch()

            else:
                print(
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

                    print(f"[{datetime.now()}] {wildcards.experiment}: Aggregation complete: {output.agg_count_txt}")

                    # --- REMOVED DIRECTORY DELETION ---
                    # intermediate_counts_dir = Path(sample_files[0]).parent
                    # ... removal logic removed ...

                except Exception as e:
                    print(
                        f"[{datetime.now()}] ERROR during count aggregation for {wildcards.experiment}: {e}"
                    )
                    with open(str(log), "w") as f:
                        f.write(f"Error: {e}")
                    raise e
        else:
            print(
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

# --- Rule to run MAGeCK RRA per contrast ---
rule run_mageck_rra_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        # Depend on the static contrast file name
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
        sif=MAGECK_SIF, # Depend on the specific SIF file
    output:
        gene_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_RRA.gene_summary.txt",
        sgrna_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_RRA.sgrna_summary.txt",
    params:
        # Format analysis options from config - keep this part
        analysis_options_str=lambda wc: format_options(
            {**{"norm-method": config.get("mageck_norm_method", "median")},
             **config.get("mageck_rra_options", {})})
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_rra_{contrast}.log",
    threads: 1
    container:
        # Directly reference the SIF path variable, converted to string
        str(MAGECK_SIF)
    run:
        # Imports needed within the run block scope
        import pandas as pd
        import shlex
        from snakemake.shell import shell # Explicit import might help, though usually implicit
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
            contrast_row = df[df['name'] == wildcards.contrast]
            if not contrast_row.empty:
                # Assuming treatment and control columns contain comma-separated strings already
                treatment_samples = contrast_row.iloc[0]['treatment']
                control_samples = contrast_row.iloc[0]['control']
            else:
                raise ValueError(f"Contrast '{wildcards.contrast}' not found in {input.contrasts_txt}")
        except Exception as e:
             # Log error and raise
            error_msg = f"Error parsing contrasts file {input.contrasts_txt} for contrast {wildcards.contrast}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f: f.write(error_msg)
            raise RuntimeError(error_msg) from e # Raise exception to stop the rule

        # Construct the command string using a relative output prefix
        # Snakemake will execute this in the mounted output directory
        command = f"""
        mageck test \\
            -k {input.count_file} \\
            -t {shlex.quote(treatment_samples)} \\
            -c {shlex.quote(control_samples)} \\
            -n {wildcards.contrast}_RRA \\
            {params.analysis_options_str} \\
            > {log} 2>&1
        """
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
        # Format analysis options from config
        analysis_options_str=lambda wc: format_options(config.get("mageck_mle_options", {})),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_mle_{experiment}.log", # Log per experiment
    threads: 1 # MAGeCK MLE can sometimes use more, but often limited by I/O
    container:
        # Directly reference the SIF path variable, converted to string
        str(MAGECK_SIF)
    shell:
        # Snakemake ensures the output directory exists and mounts it.
        # Run mageck mle using a relative output prefix.
        r"""
        mageck mle \\
            -k {input.count_file} \\
            -d {input.design_matrix} \\
            -n {wildcards.experiment}_MLE \\
            {params.analysis_options_str} \\
            > {log} 2>&1
        """


# --- Rule to run DrugZ per contrast (Conditional) ---
rule run_drugz_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        # Depend on the static contrast file name
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
        sif=DRUGZ_SIF, # Depend on the specific SIF file
    output:
        drugz_results=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_DrugZ.txt",
    params:
        # Define output path relative to container working directory
        # output_path=lambda wc, output: output.drugz_results,
        # Format analysis options from config - keep this part
        analysis_options_str=lambda wc: format_options(config.get("drugz_options", {})),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "drugz_{contrast}.log",
    threads: 1
    container:
        # Directly reference the SIF path variable, converted to string
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
            contrast_row = df[df['name'] == wildcards.contrast]
            if not contrast_row.empty:
                # Assuming treatment and control columns contain comma-separated strings already
                treatment_samples = contrast_row.iloc[0]['treatment']
                control_samples = contrast_row.iloc[0]['control']
            else:
                raise ValueError(f"Contrast '{wildcards.contrast}' not found in {input.contrasts_txt}")
        except Exception as e:
             # Log error and raise
            error_msg = f"Error parsing contrasts file {input.contrasts_txt} for contrast {wildcards.contrast}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f: f.write(error_msg)
            raise RuntimeError(error_msg) from e # Raise exception to stop the rule

        # Construct the command string using a relative output filename
        # Snakemake will execute this in the mounted output directory
        output_filename = f"{wildcards.contrast}_DrugZ.txt"
        command = f"""
        python /drugz/drugz.py \\
            --input {input.count_file} \\
            --output "{output_filename}" \\
            --control-id {shlex.quote(control_samples)} \\
            --treatment-id {shlex.quote(treatment_samples)} \\
            {params.analysis_options_str} \\
            > {log} 2>&1
        """
        # Execute the command
        shell(command)


# --- Rule to Convert MAGeCK RRA Results to CSV ---
rule convert_rra_results:
    input:
        rra_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_RRA.gene_summary.txt",
    output:
        # Place converted file directly under experiment dir
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}_gMGK.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_rra_{contrast}.log",
    run:
        if config.get("skip_rra", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping RRA result conversion for {wildcards.contrast} due to skip_rra flag."
            print(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.csv_summary)
            sys.exit(0)

        if (
            not Path(input.rra_summary).exists()
            or Path(input.rra_summary).stat().st_size == 0
        ):
            # Log error and raise (input is required for conversion)
            error_msg = f"Input RRA summary {input.rra_summary} not found or empty. Cannot convert."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: Converting MAGeCK RRA results..."
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            convert_results_to_csv(
                input_txt_path=str(input.rra_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="rra",
            )
            print(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: RRA results converted to CSV: {output.csv_summary}")
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_results_to_csv function not found/imported."
            )
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting RRA results for {wildcards.contrast}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            raise e


# --- Rule to Convert MAGeCK MLE Results to CSV (Conditional) ---
rule convert_mle_results:
    input:
        # Input is the experiment-level MLE summary
        mle_summary=OUTPUT_DIR / "{experiment}" / "MLE_analysis_results" / "{experiment}_MLE.gene_summary.txt", # Updated input path
    output:
        # Place converted file directly under experiment dir
        csv_summary=OUTPUT_DIR / "{experiment}" / "{experiment}_gMLE.csv",
    log:
        # Log is per experiment
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_mle_{experiment}.log",
    run:
        if config.get("skip_mle", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping MLE result conversion for {wildcards.experiment} due to skip_mle flag."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.csv_summary)
            sys.exit(0)

        design_matrix_path = OUTPUT_DIR / wildcards.experiment / "design_matrix.txt"
        if not design_matrix_path.exists() or design_matrix_path.stat().st_size == 0:
            # Log skip and exit cleanly (matches skip condition in run_mle rule)
            info_msg = f"Skipping MLE result conversion for {wildcards.experiment}: Design matrix was missing or empty."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.csv_summary)
            sys.exit(0)

        if (
            not Path(input.mle_summary).exists()
            or Path(input.mle_summary).stat().st_size == 0
        ):
            # Log error and raise (input is required for conversion)
            error_msg = f"Input MLE summary {input.mle_summary} not found or empty. Cannot convert."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"[{datetime.now()}] {wildcards.experiment}: Converting MAGeCK MLE results..."
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            # Ensure this function exists and handles the input/output format
            convert_results_to_csv(
                input_txt_path=str(input.mle_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="mle",
            )
            print(f"[{datetime.now()}] {wildcards.experiment}: MLE results converted to CSV: {output.csv_summary}")
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_results_to_csv function not found/imported."
            )
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting MLE results for {wildcards.experiment}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            raise e


# --- Rule to Convert DrugZ Results to CSV (Conditional) ---
rule convert_drugz_results:
    input:
        drugz_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_DrugZ.txt",
    output:
        # Place converted file directly under experiment dir
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}_gDZ.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_drugz_{contrast}.log",
    run:
        if config.get("skip_drugz", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping DrugZ result conversion for {wildcards.contrast} due to skip_drugz flag."
            print(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.csv_summary)
            sys.exit(0)

        if (
            not Path(input.drugz_summary).exists()
            or Path(input.drugz_summary).stat().st_size == 0
        ):
            # Log error and raise (input required)
            error_msg = f"Input DrugZ summary {input.drugz_summary} not found or empty. Cannot convert."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: Converting DrugZ results..."
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            convert_results_to_csv(
                input_txt_path=str(input.drugz_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="drugz",
            )
            print(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: DrugZ results converted to CSV: {output.csv_summary}")
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_results_to_csv function not found/imported."
            )
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting DrugZ results for {wildcards.contrast}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
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
        if config.get("skip_qc", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping sgRNA distribution plot for {wildcards.experiment} due to skip_qc flag."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.html_plot)
            sys.exit(0)

        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Skipping sgRNA distribution plot."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(f"[{datetime.now()}] {wildcards.experiment}: Plotting sgRNA distribution...")
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_sgRNA_distribution(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting sgRNA distribution failed: {msg}"
                print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            print(f"[{datetime.now()}] {wildcards.experiment}: sgRNA distribution plot created: {output.html_plot}")

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_sgRNA_distribution function not found/imported."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Error plotting sgRNA distribution for {wildcards.experiment}: {e}"
            )
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
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
        if config.get("skip_qc", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping gene distribution plot for {wildcards.experiment} due to skip_qc flag."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.html_plot)
            sys.exit(0)

        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Skipping gene distribution plot."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(f"[{datetime.now()}] {wildcards.experiment}: Plotting gene distribution...")
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_gene_distribution(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting gene distribution failed: {msg}"
                print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            print(f"[{datetime.now()}] {wildcards.experiment}: Gene distribution plot created: {output.html_plot}")

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_gene_distribution function not found/imported."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Error plotting gene distribution for {wildcards.experiment}: {e}"
            )
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
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
        if config.get("skip_qc", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping Gini index plot for {wildcards.experiment} due to skip_qc flag."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.html_plot)
            sys.exit(0)

        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Skipping Gini index plot."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(f"[{datetime.now()}] {wildcards.experiment}: Plotting Gini index...")
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_gini_index(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting Gini index failed: {msg}"
                print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            print(f"[{datetime.now()}] {wildcards.experiment}: Gini index plot created: {output.html_plot}")

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_gini_index function not found/imported."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error plotting Gini index for {wildcards.experiment}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            raise e


# Placeholder rule to plot ROC curve (requires controls)
rule plot_roc_curve:
    input:
        mageck_results=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "{contrast}_gMGK.csv",
        known_controls=lambda wc: BASE_DIR / wc.experiment / "known_controls.csv",
    output:
        html_plot=OUTPUT_DIR
        / "{experiment}"
        / "qc"
        / "{experiment}_{contrast}_roc.html",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "plot_roc_{contrast}.log",
    run:
        if config.get("skip_qc", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping ROC curve plot for {wildcards.experiment}/{wildcards.contrast} due to skip_qc flag."
            print(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.html_plot)
            sys.exit(0)

        if (
            not Path(input.mageck_results).exists()
            or Path(input.mageck_results).stat().st_size == 0
        ):
            # Log error and raise (required input)
            error_msg = f"Input MAGeCK results {input.mageck_results} not found or empty. Cannot plot ROC."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        if not Path(input.known_controls).exists():
            # Log skip and exit cleanly (optional input)
            warn_msg = f"Known controls file {input.known_controls} not found. Skipping ROC plot for {wildcards.experiment}/{wildcards.contrast}."
            print(f"[{datetime.now()}] WARNING {wildcards.experiment}/{wildcards.contrast}: {warn_msg}")
            with open(str(log), "w") as f:
                f.write(warn_msg)
            # touch(output.html_plot)
            sys.exit(0)

        print(
            f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: Plotting ROC curve (PLACEHOLDER)..."
        )
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            print(
                "--> Placeholder rule: plot_roc_curve() function would be called here."
            )
            # Placeholder: Just log completion, no touch needed
            # touch(output.html_plot)
            print(f"[{datetime.now()}] {wildcards.experiment}/{wildcards.contrast}: ROC curve plot logged completion (placeholder): {output.html_plot}")

        except NameError:
            error_msg = "ERROR: plot_roc_curve function not found/imported (placeholder assumes it exists)."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # No touch needed for placeholder
            # touch(output.html_plot)
        except Exception as e:
            error_msg = f"Error plotting ROC curve for {wildcards.experiment}/{wildcards.contrast}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}/{wildcards.contrast}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
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
        roc_curves=lambda wc: expand(
            OUTPUT_DIR
            / wc.experiment
            / "qc"
            / f"{wc.experiment}_{{contrast}}_roc.html",
            contrast=get_contrast_names(SimpleNamespace(experiment=wc.experiment)),
        ),
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
        if config.get("skip_qc", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping QC report generation for {wildcards.experiment} due to skip_qc flag."
            print(f"[{datetime.now()}] {wildcards.experiment}: {info_msg}")
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.report)
            sys.exit(0)

        print(f"[{datetime.now()}] {wildcards.experiment}: Generating QC report (PLACEHOLDER)...")
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

            print(
                "--> Placeholder rule: generate_qc_report() function would be called here with:"
            )
            print(input_files_dict)
            # Placeholder: Log completion, no touch needed
            # touch(output.report)
            print(f"[{datetime.now()}] {wildcards.experiment}: QC report logged completion (placeholder): {output.report}")

        except NameError:
            error_msg = "ERROR: generate_qc_report function not found/imported (placeholder assumes it exists)."
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # No touch needed
            # touch(output.report)
        except Exception as e:
            error_msg = f"Error generating QC report for {wildcards.experiment}: {e}"
            print(f"[{datetime.now()}] ERROR {wildcards.experiment}: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.report)
            raise e


# --- Rule All: Define Final Output Files ---


def get_final_outputs(wildcards, checkpoints):
    """Dynamically collects all expected final output files based on config and targets."""
    final_files = []
    print("--- Entering get_final_outputs ---") # DEBUG
    print(f"Target experiments: {TARGET_EXPERIMENTS}") # DEBUG

    # Use wildcards.experiment if rule all defines experiment, otherwise iterate
    # Assuming rule all doesn't have {experiment} wildcard for simplicity now
    # We iterate over TARGET_EXPERIMENTS determined earlier
    for experiment in TARGET_EXPERIMENTS:
        print(f"Processing experiment: {experiment}") # DEBUG
        validation_info = get_validation_info(experiment)
        print(f"Validation info for {experiment}: {validation_info}") # DEBUG

        if validation_info["status"] == "failed":
            print(
                f"Skipping final output collection for failed experiment: {experiment}"
            )
            continue

        # Use the modified get_contrast_names function with checkpoints object
        # Pass SimpleNamespace for wildcards as before if needed, or adjust based on how rule all calls this
        print(f"Calling get_contrast_names for {experiment}...") # DEBUG
        contrasts = get_contrast_names(SimpleNamespace(experiment=experiment), checkpoints)
        print(f"Contrasts for {experiment}: {contrasts}") # DEBUG

        if (
            not isinstance(contrasts, list) or (contrasts and "error" in contrasts[0]) # Handle error case and empty list
        ):
            print(
                f"Skipping contrast-specific outputs for {experiment} due to contrast parsing issues."
            )
            contrasts = [] # Ensure contrasts is an empty list if issues occurred

        # --- 1. Analysis Results (per contrast / per experiment) ---
        print(f"Checking analysis results for {experiment}...") # DEBUG
        for contrast in contrasts:
            print(f"  Checking contrast: {contrast}") # DEBUG
            # Add RRA results if not skipped
            if not config.get("skip_rra", False):
                rra_file = OUTPUT_DIR / experiment / f"{contrast}_gMGK.csv"
                print(f"    Adding RRA target: {rra_file}") # DEBUG
                final_files.append(rra_file)
            else:
                print("    Skipping RRA (skip_rra=True)") # DEBUG

            # Add DrugZ results if not skipped
            if not config.get("skip_drugz", False):
                dz_file = OUTPUT_DIR / experiment / f"{contrast}_gDZ.csv"
                print(f"    Adding DrugZ target: {dz_file}") # DEBUG
                final_files.append(dz_file)
            else:
                 print("    Skipping DrugZ (skip_drugz=True)") # DEBUG

        # Add MLE results (per experiment) if not skipped AND design matrix exists
        print(f"Checking MLE results for {experiment}...") # DEBUG
        # Check for *converted* design matrix in output dir
        design_matrix_path = OUTPUT_DIR / experiment / "design_matrix.txt"
        # Check if the *rule output* design matrix exists and is non-empty
        design_output_exists = design_matrix_path.exists() # and design_matrix_path.stat().st_size > 0 # Size check might be too strict initially
        print(f"  MLE Check: skip_mle={config.get('skip_mle', False)}, design_matrix_exists={design_output_exists} ({design_matrix_path})") # DEBUG

        if not config.get("skip_mle", False) and design_output_exists:
             # Add NATIVE MLE outputs (now experiment-level)
            mle_gene_native = OUTPUT_DIR / experiment / "analysis_results" / f"{experiment}_MLE.gene_summary.txt"
            mle_sgrna_native = OUTPUT_DIR / experiment / "analysis_results" / f"{experiment}_MLE.sgrna_summary.txt"
            mle_beta_native = OUTPUT_DIR / experiment / "analysis_results" / f"{experiment}_MLE.beta_coefficients.txt"
            # Add CONVERTED MLE CSV (at experiment level)
            mle_csv = OUTPUT_DIR / experiment / f"{experiment}_gMLE.csv"
            print(f"    Adding MLE targets: {mle_gene_native}, {mle_sgrna_native}, {mle_beta_native}, {mle_csv}") # DEBUG
            final_files.append(mle_gene_native)
            final_files.append(mle_sgrna_native)
            final_files.append(mle_beta_native)
            final_files.append(mle_csv)
        else:
            print(f"    Skipping MLE targets.") # DEBUG

        # --- 2. QC Files (Conditional) ---
        print(f"Checking QC files for {experiment}...") # DEBUG
        if not config.get("skip_qc", False):
            # Per-experiment QC plots (ensure rules still exist)
            if "plot_sgrna_distribution" in globals(): # Check if rule exists
                sgrna_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_sgrna_distribution.html"
                print(f"    Adding QC target: {sgrna_plot}") # DEBUG
                final_files.append(sgrna_plot)
            if "plot_gene_distribution" in globals():
                gene_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_gene_distribution.html"
                print(f"    Adding QC target: {gene_plot}") # DEBUG
                final_files.append(gene_plot)
            if "plot_gini_index" in globals():
                gini_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_gini_index.html"
                print(f"    Adding QC target: {gini_plot}") # DEBUG
                final_files.append(gini_plot)
            if "generate_qc_report" in globals():
                qc_report = OUTPUT_DIR / experiment / "qc" / f"{experiment}_qc_report.html"
                print(f"    Adding QC target: {qc_report}") # DEBUG
                final_files.append(qc_report)

            # Per-contrast QC plots (ROC needs controls check)
            if "plot_roc_curve" in globals():
                known_controls_path = BASE_DIR / experiment / "known_controls.csv"
                print(f"  ROC Check: known_controls_exists={known_controls_path.exists()} ({known_controls_path})") # DEBUG
                if known_controls_path.exists(): # Only add ROC if controls file exists
                    for contrast in contrasts:
                        roc_plot = OUTPUT_DIR / experiment / "qc" / f"{experiment}_{contrast}_roc.html"
                        print(f"    Adding QC target: {roc_plot}") # DEBUG
                        final_files.append(roc_plot)
                else:
                    print("    Skipping ROC plots (no known_controls.csv)") # DEBUG
            else:
                 print("    Skipping ROC plots (rule not defined)") # DEBUG

            # Per-sample FastQC reports (if FASTQ input)
            print(f"  FastQC Check: data_type={validation_info.get('data_type')}") # DEBUG
            if validation_info.get("data_type") == "fastq" and "run_fastqc_per_sample" in globals():
                fastq_basenames = get_fastq_basenames(experiment)
                print(f"    FastQC basenames: {fastq_basenames}") # DEBUG
                for sample in fastq_basenames:
                    fq_html = OUTPUT_DIR / experiment / "qc" / f"{sample}_fastqc.html"
                    fq_zip = OUTPUT_DIR / experiment / "qc" / f"{sample}_fastqc.zip"
                    print(f"    Adding QC targets: {fq_html}, {fq_zip}") # DEBUG
                    final_files.append(fq_html)
                    final_files.append(fq_zip)
            else:
                print("    Skipping FastQC reports (not FASTQ data or rule undefined)") # DEBUG
        else:
            print("  Skipping ALL QC (skip_qc=True)") # DEBUG

    # Convert Path objects to strings for Snakemake input list
    final_files_str = [str(f) for f in final_files]
    print(f"--- Exiting get_final_outputs with {len(final_files_str)} targets: {final_files_str} ---") # DEBUG
    return final_files_str


# --- Final Target Rule ---
# Remove the temporary get_first_experiment_qc_targets function
# def get_first_experiment_qc_targets(): ...


rule all:
    input:
        # Use a lambda function to pass checkpoints object to get_final_outputs
        lambda wildcards, checkpoints: get_final_outputs(wildcards, checkpoints),
    output:
        # Single flag file indicating completion of requested targets
        OUTPUT_DIR / "pipeline_complete.flag",
    run:
        print(f"[{datetime.now()}] CRISPR Analysis Pipeline Workflow Completed for targets:")
        # Input might be complex now, iterate carefully or just log completion
        # for f in input:
        #     print(f"- {f}")
        print(f"Input files generated (list might be long). Check output directory.")
        print(f"[{datetime.now()}] Completion flag: {output}")
        # Create the flag file
        Path(output[0]).touch()
