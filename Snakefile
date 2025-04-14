# flake8: noqa
# pylint: disable=all
# type: ignore

"""
Snakefile for CRISPR Analysis Pipeline (Refactored)
"""

import os
import glob
from pathlib import Path
import sys  # Added for sys.exit
import pandas as pd  # Needed for contrast parsing
import shutil  # Added for shutil.copy
from types import SimpleNamespace  # Needed for mocking wildcards

# --- Core Function Imports ---
# Assuming core modules are importable relative to Snakefile location or via PYTHONPATH
# Adjust paths as needed based on project structure
try:
    # Use absolute import path assuming Snakefile is at root
    from core.validation import validate_experiment_structure
    from core.file_handling import convert_file_to_tab_delimited, parse_contrasts
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
config.setdefault("skip_mle", False)
config.setdefault("skip_drugz", False)
config.setdefault("use_apptainer", True)  # Default to Apptainer
# Add defaults for other analysis options (e.g., mageck norm_method) if needed
config.setdefault("mageck_norm_method", "median")

# --- Helper Variables ---
BASE_DIR = Path(config["base_dir"]).resolve()
OUTPUT_DIR = Path(config["output_dir"]).resolve()
print(f"Base Input Directory: {BASE_DIR}")
print(f"Base Output Directory: {OUTPUT_DIR}")

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
# This rule acts as an implicit validation trigger
rule convert_contrasts:
    input:
        # Define potential input CSV path
        csv=lambda wc: str(BASE_DIR / wc.experiment / "contrasts.csv"),
    output:
        # Place converted file directly in experiment output dir
        txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_contrasts.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)
        if validation_info["status"] == "failed":
            error_msg = f"Skipping contrast conversion for {wildcards.experiment} due to validation failure: {validation_info.get('error')}"
            print(error_msg)
            with open(log[0], "w") as f:
                f.write(error_msg)
        else:
            print(f"Converting contrasts for {wildcards.experiment}")
            # Output dir is the experiment dir
            output_dir_path = Path(output.txt).parent
            output_dir_path.mkdir(parents=True, exist_ok=True)
            try:
                input_csv_path = validation_info["contrasts_path"]
                if Path(input_csv_path).suffix.lower() == ".csv":
                    # Pass the specific output path to the function
                    convert_file_to_tab_delimited(
                        file_path=input_csv_path,
                        output_path=str(output.txt),  # Specify exact output path
                    )
                else:
                    print(
                        f"""Input contrast file {input_csv_path} is already .txt, copying."""
                    )
                    shutil.copy(input_csv_path, output.txt)
            except Exception as e:
                print(f"Error converting contrasts for {wildcards.experiment}: {e}")
                with open(log[0], "w") as f:
                    f.write(f"Error: {e}")
                raise e  # Fail the job


# --- Rule to Convert Design Matrix CSV to TXT (Conditional) ---
rule convert_design_matrix:
    input:
        # Define potential input CSV path
        csv=lambda wc: str(BASE_DIR / wc.experiment / "design_matrix.csv"),
    output:
        # Place converted file directly in experiment output dir
        txt=OUTPUT_DIR / "{experiment}" / "design_matrix.txt",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_design_matrix.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)
        if validation_info["status"] == "failed":
            # Validation failed for the experiment, just log and keep placeholder
            print(
                f"Skipping design matrix conversion for {wildcards.experiment} due to validation failure."
            )
        elif validation_info.get("design_matrix_path"):
            # Design matrix exists, proceed with conversion
            input_path = validation_info["design_matrix_path"]
            output_dir_path = Path(output.txt).parent
            output_dir_path.mkdir(parents=True, exist_ok=True)
            if Path(input_path).suffix.lower() == ".csv":
                print(f"Converting design matrix for {wildcards.experiment}")
                try:
                    convert_file_to_tab_delimited(
                        file_path=input_path,
                        output_path=str(output.txt),  # Specify exact output path
                    )
                except Exception as e:
                    print(
                        f"Error converting design matrix for {wildcards.experiment}: {e}"
                    )
                    with open(log[0], "w") as f:
                        f.write(f"Error: {e}")
                    raise e
            else:
                print(
                    f"Input design matrix file {input_path} is already .txt, copying."
                )
                shutil.copy(input_path, output.txt)
        else:
            print(
                f"Design matrix not found for {wildcards.experiment}, skipping conversion."
            )


# --- Rule to Convert Input Read Count CSV to TSV ---
rule convert_read_count_input:
    input:
        rc_csv=lambda wc: get_validation_info(wc.experiment).get("rc_path"),
    output:
        # Output is still defined/touched by Snakemake initially
        count_txt=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_read_count.log",
    run:
        validation_info = get_validation_info(wildcards.experiment)

        if (
            validation_info["status"] == "valid"
            and validation_info.get("data_type") == "rc"
        ):
            input_rc_path = validation_info.get("rc_path")
            if input_rc_path and Path(input_rc_path).exists():
                print(
                    f"Converting input read count file for {wildcards.experiment}: {input_rc_path} to {output.count_txt}"
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
                        print(f"Renaming {created_file_path} to {target_file_path}")
                        created_file_path.rename(target_file_path)
                    elif not target_file_path.exists():
                        if not created_file_path.exists():
                            # Let the error propagate if file not created
                            raise FileNotFoundError(
                                f"make_count_table did not produce expected file: {created_file_path} or {target_file_path}"
                            )
                    print(f"Read count conversion complete: {output.count_txt}")
                except Exception as e:
                    # Log error and raise
                    error_msg = f"Error converting input read count for {wildcards.experiment}: {e}"
                    print(f"ERROR: {error_msg}")
                    with open(str(log), "w") as f:
                        f.write(error_msg)
                    raise e
            else:
                # Log warning and exit cleanly (no output expected)
                warn_msg = f"Warning: data_type is 'rc' but rc_path not found/valid for {wildcards.experiment}. Skipping conversion."
                print(warn_msg)
                with open(str(log), "w") as f:
                    f.write(warn_msg)
                # No touch() needed here, exit cleanly
                sys.exit(0)
        else:
            # Log info and exit cleanly (no output expected)
            info_msg = f"Skipping read count conversion for {wildcards.experiment} (data type: {validation_info.get('data_type')}, status: {validation_info.get('status')})."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # No touch() needed here, exit cleanly
            sys.exit(0)


# --- Function to parse contrasts from the *converted* TXT file ---
ContrastCache = {}


def get_contrast_names(wildcards):
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
            return ContrastCache[cache_key]

    # Path to the *output* of the conversion rule
    converted_contrast_path = OUTPUT_DIR / experiment / "contrasts.txt"

    # We need the conversion rule to run first.
    # This function might be called by 'rule all' before the file exists.
    # A better approach might be a checkpoint for parsing.
    # For now, rely on rule dependency and check existence.
    validation_info = get_validation_info(experiment)
    if validation_info["status"] == "failed":
        ContrastCache[cache_key] = ["validation_failed"]
        return ["validation_failed"]

    if not converted_contrast_path.exists():
        print(
            f"Warning: Waiting for converted contrasts file: {converted_contrast_path}"
        )
        # Snakemake should handle this dependency via rule all input
        # If called elsewhere prematurely, this might be an issue.
        ContrastCache[cache_key] = ["conversion_pending"]
        return ["conversion_pending"]
    try:
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
    """Find the R1 or single-end FASTQ for a given sample basename."""
    validation_info = get_validation_info(wildcards.experiment)
    if (
        validation_info["status"] != "valid"
        or validation_info.get("data_type") != "fastq"
    ):
        return None  # Or raise error

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
    return None  # Not found


# --- Function to get the R2 FASTQ file for a sample basename ---
def get_fastq_r2_for_sample(wildcards):
    """Find the R2 FASTQ for a given sample basename, if it exists."""
    validation_info = get_validation_info(wildcards.experiment)
    if (
        validation_info["status"] != "valid"
        or validation_info.get("data_type") != "fastq"
    ):
        return None

    fastq_files = validation_info.get("fastq_files", [])
    # Look for R2 specifically
    for f in fastq_files:
        if Path(f).name.startswith(wildcards.sample) and "_R2" in Path(f).name:
            return f
    return None  # R2 not found


# --- Rule to run FastQC on individual FASTQ files ---
rule run_fastqc_per_sample:
    input:
        fastq=get_fastq_for_sample,
    output:
        html=OUTPUT_DIR / "{experiment}" / "qc" / "{sample}_fastqc.html",
        zip=OUTPUT_DIR / "{experiment}" / "qc" / "{sample}_fastqc.zip",
    params:
        output_dir=lambda wc, output: str(Path(output.html).parent),
        use_apptainer=config["use_apptainer"],
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "fastqc_{sample}.log",
    threads: 1
    run:
        if not input.fastq:
            # Raise error for missing required input
            error_msg = f"Input FASTQ file could not be determined for sample '{wildcards.sample}' in experiment '{wildcards.experiment}'. Check validation."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise FileNotFoundError(error_msg)

        if config["skip_qc"]:
            # Log skip and exit cleanly, no touch() needed
            info_msg = f"Skipping FastQC for {wildcards.sample} due to skip_qc flag."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.html)
            # touch(output.zip)
            sys.exit(0)
        else:
            print(
                f"Running FastQC for {wildcards.experiment} - {wildcards.sample} ({input.fastq})"
            )
            Path(params.output_dir).mkdir(parents=True, exist_ok=True)
            try:
                success, out_path_or_msg = run_fastqc(
                    fastq_path=str(input.fastq),
                    output_dir=params.output_dir,
                    use_apptainer=params.use_apptainer,
                )
                if not success:
                    # Log error and raise
                    error_msg = (
                        f"run_fastqc failed for {input.fastq}: {out_path_or_msg}"
                    )
                    print(f"ERROR: {error_msg}")
                    with open(str(log), "w") as f:
                        f.write(error_msg)
                    raise RuntimeError(error_msg)
                print(f"FastQC finished successfully for {wildcards.sample}")
            except Exception as e:
                # Log error and raise
                error_msg = f"Error running FastQC for {wildcards.sample}: {e}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(f"Error running FastQC on {input.fastq}:\n{e}")
                raise e


# --- Rule to Aggregate Sample Counts (from FASTQ processing) ---
rule aggregate_counts:
    input:
        sample_counts=lambda wc: expand(
            OUTPUT_DIR / "{experiment}" / "counts" / "{sample}.count.txt",
            experiment=wc.experiment,
            sample=get_fastq_basenames(wc.experiment),
        ),
    output:
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
                    f"Warning: No sample count files found to aggregate for {wildcards.experiment}. Creating empty output."
                )
            else:
                print(
                    f"Aggregating {len(sample_files)} sample counts for {wildcards.experiment}"
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

                    print(f"Aggregation complete: {output.agg_count_txt}")

                    # Delete intermediate counts directory on success
                    intermediate_counts_dir = Path(sample_files[0]).parent
                    if (
                        intermediate_counts_dir.exists()
                        and intermediate_counts_dir.name == "counts"
                    ):
                        print(
                            f"Removing intermediate counts directory: {intermediate_counts_dir}"
                        )
                        shutil.rmtree(intermediate_counts_dir)
                    else:
                        print(
                            f"Warning: Could not reliably determine intermediate counts directory for removal ({intermediate_counts_dir}). Skipping removal."
                        )

                except Exception as e:
                    print(
                        f"Error during count aggregation or cleanup for {wildcards.experiment}: {e}"
                    )
                    with open(str(log), "w") as f:
                        f.write(f"Error: {e}")
                    raise e
        else:
            print(
                f"Skipping count aggregation for {wildcards.experiment} (data type: {validation_info.get('data_type')}, status: {validation_info.get('status')})."
            )


# --- Rule to run MAGeCK count per sample ---
rule run_mageck_count_per_sample:
    input:
        # Gets R1 or single-end
        r1=get_fastq_for_sample,
        # Gets R2 if present, otherwise None
        r2=get_fastq_r2_for_sample,
        # Library path comes from validation info, accessed via params
    output:
        # Place results in a dedicated intermediate counts directory
        count=OUTPUT_DIR / "{experiment}" / "counts" / "{sample}.count.txt",
        summary=OUTPUT_DIR / "{experiment}" / "counts" / "{sample}.countsummary.txt",
    params:
        # Get library path from validation cache
        library=lambda wc: get_validation_info(wc.experiment).get("library_path"),
        output_prefix=lambda wc, output: str(Path(output.count).parent / wc.sample),
        sample_name="{sample}",
        use_apptainer=config["use_apptainer"],
        # Pass other MAGeCK count options from config if needed
        # count_options = config.get("mageck_count_options", {})
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_count_{sample}.log",
    # Adjust if MAGeCK count can use more threads
    threads: 1
    run:
        # Ensure library was found during validation
        if not params.library:
            raise ValueError(
                f"Library path not found during validation for experiment '{wildcards.experiment}'"
            )
        if not Path(params.library).exists():
            raise FileNotFoundError(
                f"Library file {params.library} not found for experiment '{wildcards.experiment}'"
            )

            # Ensure R1 input was found
        if not input.r1:
            raise FileNotFoundError(
                f"Input R1 FASTQ file could not be determined for sample '{wildcards.sample}' in experiment '{wildcards.experiment}'."
            )

        print(f"Running MAGeCK count for {wildcards.experiment} - {wildcards.sample}")
        # Ensure output directory exists
        Path(output.count).parent.mkdir(parents=True, exist_ok=True)

        try:
            # Call the core analysis step function
            success, out_path_or_msg = run_mageck_count(
                r1_fastq=str(input.r1),
                library_path=str(params.library),
                output_prefix=str(params.output_prefix),
                r2_fastq=str(input.r2) if input.r2 else None,
                sample_name=params.sample_name,
                # count_options=params.count_options, # Pass if defined
                use_apptainer=params.use_apptainer,
            )
            if not success:
                raise RuntimeError(
                    f"run_mageck_count failed for {input.r1}: {out_path_or_msg}"
                )
            print(f"MAGeCK count finished successfully for {wildcards.sample}")
        except Exception as e:
            print(f"Error running MAGeCK count for {wildcards.sample}: {e}")
            with open(str(log), "w") as f:
                f.write(f"Error running MAGeCK count on {input.r1}:\n{e}")
            raise e


# --- Rule to run MAGeCK RRA per contrast ---
rule run_mageck_rra_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
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
        output_prefix=lambda wc, output: str(
            Path(output.gene_summary).parent / f"{wc.contrast}_RRA"
        ),
        use_apptainer=config["use_apptainer"],
        analysis_options={"norm-method": config.get("mageck_norm_method", "median")},
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_rra_{contrast}.log",
    threads: 1
    run:
        config.setdefault("skip_rra", False)
        if config["skip_rra"]:
            # Log skip and exit cleanly
            info_msg = f"Skipping MAGeCK RRA for {wildcards.experiment} - Contrast: {wildcards.contrast} due to skip_rra flag."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            sys.exit(0)

        experiment_contrasts = ContrastCache.get(wildcards.experiment)
        if not isinstance(experiment_contrasts, list):
            # Raise error
            error_msg = f"Contrast information not available or invalid for experiment '{wildcards.experiment}'."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise ValueError(error_msg)

        contrast_info = next(
            (c for c in experiment_contrasts if c.get("name") == wildcards.contrast),
            None,
        )
        if not contrast_info:
            # Raise error
            error_msg = f"Could not find contrast details for '{wildcards.contrast}' in experiment '{wildcards.experiment}'."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise ValueError(error_msg)

        print(
            f"Running MAGeCK RRA for {wildcards.experiment} - Contrast: {wildcards.contrast}"
        )

        # Check required input count file
        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Cannot run RRA for {wildcards.contrast}."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            # sys.exit(0) # Don't exit cleanly, this is an error
            raise FileNotFoundError(error_msg)

        Path(output.gene_summary).parent.mkdir(parents=True, exist_ok=True)

        try:
            success, out_path_or_msg = run_mageck_rra(
                count_path=str(input.count_file),
                contrast_info=contrast_info,
                output_prefix=params.output_prefix,
                analysis_options=params.analysis_options,
                use_apptainer=params.use_apptainer,
            )

            if not success:
                # Log error and raise
                error_msg = f"run_mageck_rra failed for contrast {wildcards.contrast}: {out_path_or_msg}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                # touch(output.gene_summary)
                # touch(output.sgrna_summary)
                raise RuntimeError(error_msg)

            expected_output = Path(out_path_or_msg)
            if str(expected_output.resolve()) != str(
                Path(output.gene_summary).resolve()
            ):
                warn_msg = f"Warning: run_mageck_rra reported success path '{out_path_or_msg}' which doesn't match expected output '{output.gene_summary}'"
                print(warn_msg)

            print(
                f"MAGeCK RRA finished successfully for {wildcards.experiment} - {wildcards.contrast}"
            )

            sgrna_path_expected = Path(output.sgrna_summary)
            sgrna_path_actual = Path(params.output_prefix + ".sgrna_summary.txt")
            if not sgrna_path_actual.exists() or str(
                sgrna_path_actual.resolve()
            ) != str(sgrna_path_expected.resolve()):
                # Just log warning about secondary output
                print(
                    f"Warning: Expected sgrna summary file {sgrna_path_expected} not found or path mismatch (actual: {sgrna_path_actual})."
                )
                # No touch() needed

        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Unexpected error running MAGeCK RRA for {wildcards.contrast}: {e}"
            )
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            raise e


# --- Rule to run MAGeCK MLE per contrast (Conditional) ---
rule run_mageck_mle_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        design_matrix=OUTPUT_DIR / "{experiment}" / "design_matrix.txt",
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
    output:
        gene_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_MLE.gene_summary.txt",
        sgrna_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_MLE.sgrna_summary.txt",
        beta_coeff=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_MLE.beta_coefficients.txt",
    params:
        output_prefix=lambda wc, output: str(
            Path(output.gene_summary).parent / f"{wc.contrast}_MLE"
        ),
        use_apptainer=config["use_apptainer"],
        analysis_options=config.get("mageck_mle_options", {}),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "mageck_mle_{contrast}.log",
    threads: 1
    run:
        if config["skip_mle"]:
            # Log skip and exit cleanly
            info_msg = f"Skipping MAGeCK MLE for {wildcards.experiment} - Contrast: {wildcards.contrast} due to skip_mle flag."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            # touch(output.beta_coeff)
            sys.exit(0)

        design_matrix_path = Path(input.design_matrix)
        if not design_matrix_path.exists() or design_matrix_path.stat().st_size == 0:
            # Log skip for missing optional input and exit cleanly
            info_msg = f"Skipping MAGeCK MLE for {wildcards.experiment} - Contrast: {wildcards.contrast}: Design matrix '{design_matrix_path}' not found or is empty."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            # touch(output.beta_coeff)
            sys.exit(0)

        experiment_contrasts = ContrastCache.get(wildcards.experiment)
        if not isinstance(experiment_contrasts, list):
            # Raise error
            error_msg = f"Contrast information not available or invalid for experiment '{wildcards.experiment}' needed for MLE rule setup."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise ValueError(error_msg)

        contrast_info = next(
            (c for c in experiment_contrasts if c.get("name") == wildcards.contrast),
            None,
        )
        if not contrast_info:
            print(
                f"Warning: Could not find exact contrast details for '{wildcards.contrast}' in parsed contrasts. Using basic name for MLE."
            )
            contrast_info = {"name": wildcards.contrast}

            # Check required input count file
        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Cannot run MLE for {wildcards.contrast}."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            # touch(output.beta_coeff)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"Running MAGeCK MLE for {wildcards.experiment} - Contrast: {wildcards.contrast}"
        )
        Path(output.gene_summary).parent.mkdir(parents=True, exist_ok=True)

        try:
            success, out_path_or_msg = run_mageck_mle(
                count_path=str(input.count_file),
                design_matrix_path=str(design_matrix_path),
                contrast_info=contrast_info,
                output_prefix=params.output_prefix,
                analysis_options=params.analysis_options,
                use_apptainer=params.use_apptainer,
            )

            if not success:
                # Log error and raise
                error_msg = f"run_mageck_mle failed for contrast {wildcards.contrast}: {out_path_or_msg}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                # touch(output.gene_summary)
                # touch(output.sgrna_summary)
                # touch(output.beta_coeff)
                raise RuntimeError(error_msg)

            expected_output = Path(out_path_or_msg)
            if str(expected_output.resolve()) != str(
                Path(output.gene_summary).resolve()
            ):
                warn_msg = f"Warning: run_mageck_mle reported success path '{out_path_or_msg}' which doesn't match expected output '{output.gene_summary}'"
                print(warn_msg)

            print(
                f"MAGeCK MLE finished successfully for {wildcards.experiment} - {wildcards.contrast}"
            )

            sgrna_path_expected = Path(output.sgrna_summary)
            sgrna_path_actual = Path(params.output_prefix + ".sgrna_summary.txt")
            if not sgrna_path_actual.exists() or str(
                sgrna_path_actual.resolve()
            ) != str(sgrna_path_expected.resolve()):
                print(
                    f"Warning: Expected sgrna summary file {sgrna_path_expected} not found or path mismatch (actual: {sgrna_path_actual})."
                )
                # No touch()

            beta_path_expected = Path(output.beta_coeff)
            beta_path_actual = Path(params.output_prefix + ".beta_coefficients.txt")
            if not beta_path_actual.exists() or str(beta_path_actual.resolve()) != str(
                beta_path_expected.resolve()
            ):
                print(
                    f"Warning: Expected beta coefficients file {beta_path_expected} not found or path mismatch (actual: {beta_path_actual})."
                )
                # No touch()

        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Unexpected error running MAGeCK MLE for {wildcards.contrast}: {e}"
            )
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.gene_summary)
            # touch(output.sgrna_summary)
            # touch(output.beta_coeff)
            raise e


# --- Rule to run DrugZ per contrast (Conditional) ---
rule run_drugz_per_contrast:
    input:
        count_file=OUTPUT_DIR / "{experiment}" / "{experiment}.count.txt",
        contrasts_txt=OUTPUT_DIR / "{experiment}" / "contrasts.txt",
    output:
        drugz_results=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_DrugZ.txt",
    params:
        output_prefix=lambda wc, output: str(
            Path(output.drugz_results).parent / f"{wc.contrast}_DrugZ"
        ),
        use_apptainer=config["use_apptainer"],
        analysis_options=config.get("drugz_options", {}),
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "drugz_{contrast}.log",
    threads: 1
    run:
        if config["skip_drugz"]:
            # Log skip and exit cleanly
            info_msg = f"Skipping DrugZ for {wildcards.experiment} - Contrast: {wildcards.contrast} due to skip_drugz flag."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.drugz_results)
            sys.exit(0)

        experiment_contrasts = ContrastCache.get(wildcards.experiment)
        if not isinstance(experiment_contrasts, list):
            # Raise error
            error_msg = f"Contrast information not available or invalid for experiment '{wildcards.experiment}'."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise ValueError(error_msg)

        contrast_info = next(
            (c for c in experiment_contrasts if c.get("name") == wildcards.contrast),
            None,
        )
        if not contrast_info:
            # Raise error
            error_msg = f"Could not find contrast details for '{wildcards.contrast}' in experiment '{wildcards.experiment}'."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise ValueError(error_msg)

            # Check required input count file
        if (
            not Path(input.count_file).exists()
            or Path(input.count_file).stat().st_size == 0
        ):
            # Log error and raise
            error_msg = f"Input count file {input.count_file} is missing or empty. Cannot run DrugZ for {wildcards.contrast}."
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.drugz_results)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"Running DrugZ for {wildcards.experiment} - Contrast: {wildcards.contrast}"
        )
        Path(output.drugz_results).parent.mkdir(parents=True, exist_ok=True)

        try:
            success, out_path_or_msg = run_drugz(
                count_path=str(input.count_file),
                contrast_info=contrast_info,
                output_prefix=params.output_prefix,
                analysis_options=params.analysis_options,
                use_apptainer=params.use_apptainer,
            )

            if not success:
                # Log error and raise
                error_msg = f"run_drugz failed for contrast {wildcards.contrast}: {out_path_or_msg}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                # touch(output.drugz_results)
                raise RuntimeError(error_msg)

            expected_output = Path(out_path_or_msg)
            if str(expected_output.resolve()) != str(
                Path(output.drugz_results).resolve()
            ):
                warn_msg = f"Warning: run_drugz reported success path '{out_path_or_msg}' which doesn't match expected output '{output.drugz_results}'"
                print(warn_msg)

            print(
                f"DrugZ finished successfully for {wildcards.experiment} - {wildcards.contrast}"
            )

        except Exception as e:
            # Log error and raise
            error_msg = f"Unexpected error running DrugZ for {wildcards.contrast}: {e}"
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.drugz_results)
            raise e


# --- Rule to Convert MAGeCK RRA Results to CSV ---
rule convert_rra_results:
    input:
        rra_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_RRA.gene_summary.txt",
    output:
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}" / "{contrast}_gMGK.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_rra_{contrast}.log",
    run:
        if config.get("skip_rra", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping RRA result conversion for {wildcards.contrast} due to skip_rra flag."
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"Converting MAGeCK RRA results for {wildcards.experiment} - Contrast: {wildcards.contrast}"
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            convert_mageck_summary_to_csv(
                input_txt_path=str(input.rra_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="rra",
            )
            print(f"RRA results converted to CSV: {output.csv_summary}")
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_mageck_summary_to_csv function not found/imported."
            )
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting RRA results for {wildcards.contrast}: {e}"
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            raise e


# --- Rule to Convert MAGeCK MLE Results to CSV (Conditional) ---
rule convert_mle_results:
    input:
        mle_summary=OUTPUT_DIR
        / "{experiment}"
        / "{contrast}"
        / "analysis_results"
        / "{contrast}_MLE.gene_summary.txt",
    output:
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}" / "{contrast}_gMLE.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_mle_{contrast}.log",
    run:
        if config.get("skip_mle", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping MLE result conversion for {wildcards.contrast} due to skip_mle flag."
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.csv_summary)
            sys.exit(0)

        design_matrix_path = OUTPUT_DIR / wildcards.experiment / "design_matrix.txt"
        if not design_matrix_path.exists() or design_matrix_path.stat().st_size == 0:
            # Log skip and exit cleanly (matches skip condition in run_mle rule)
            info_msg = f"Skipping MLE result conversion for {wildcards.contrast}: Design matrix was missing or empty."
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"Converting MAGeCK MLE results for {wildcards.experiment} - Contrast: {wildcards.contrast}"
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            convert_mageck_summary_to_csv(
                input_txt_path=str(input.mle_summary),
                output_csv_path=str(output.csv_summary),
                analysis_type="mle",
            )
            print(f"MLE results converted to CSV: {output.csv_summary}")
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_mageck_summary_to_csv function not found/imported."
            )
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting MLE results for {wildcards.contrast}: {e}"
            print(f"ERROR: {error_msg}")
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
        csv_summary=OUTPUT_DIR / "{experiment}" / "{contrast}" / "{contrast}_gDZ.csv",
    log:
        OUTPUT_DIR / "{experiment}" / "logs" / "convert_drugz_{contrast}.log",
    run:
        if config.get("skip_drugz", False):
            # Log skip and exit cleanly
            info_msg = f"Skipping DrugZ result conversion for {wildcards.contrast} due to skip_drugz flag."
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.csv_summary)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(
            f"Converting DrugZ results for {wildcards.experiment} - Contrast: {wildcards.contrast}"
        )
        try:
            Path(output.csv_summary).parent.mkdir(parents=True, exist_ok=True)
            convert_drugz_summary_to_csv(
                input_txt_path=str(input.drugz_summary),
                output_csv_path=str(output.csv_summary),
            )
            print(f"DrugZ results converted to CSV: {output.csv_summary}")
        except NameError:
            # Log error and raise
            error_msg = (
                "ERROR: convert_drugz_summary_to_csv function not found/imported."
            )
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error converting DrugZ results for {wildcards.contrast}: {e}"
            print(f"ERROR: {error_msg}")
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
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(f"Plotting sgRNA distribution for {wildcards.experiment}")
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_sgRNA_distribution(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting sgRNA distribution failed: {msg}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            print(f"sgRNA distribution plot created: {output.html_plot}")

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_sgRNA_distribution function not found/imported."
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Error plotting sgRNA distribution for {wildcards.experiment}: {e}"
            )
            print(f"ERROR: {error_msg}")
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
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(f"Plotting gene distribution for {wildcards.experiment}")
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_gene_distribution(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting gene distribution failed: {msg}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            print(f"Gene distribution plot created: {output.html_plot}")

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_gene_distribution function not found/imported."
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = (
                f"Error plotting gene distribution for {wildcards.experiment}: {e}"
            )
            print(f"ERROR: {error_msg}")
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
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        print(f"Plotting Gini index for {wildcards.experiment}")
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            success, msg = plot_gini_index(
                count_path=str(input.count_file), output_html_path=str(output.html_plot)
            )

            if not success:
                # Log error and raise
                error_msg = f"Plotting Gini index failed: {msg}"
                print(f"ERROR: {error_msg}")
                with open(str(log), "w") as f:
                    f.write(error_msg)
                raise RuntimeError(error_msg)

            print(f"Gini index plot created: {output.html_plot}")

        except NameError:
            # Log error and raise
            error_msg = "ERROR: plot_gini_index function not found/imported."
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            raise
        except Exception as e:
            # Log error and raise
            error_msg = f"Error plotting Gini index for {wildcards.experiment}: {e}"
            print(f"ERROR: {error_msg}")
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
            print(info_msg)
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
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.html_plot)
            # sys.exit(0)
            raise FileNotFoundError(error_msg)

        if not Path(input.known_controls).exists():
            # Log skip and exit cleanly (optional input)
            warn_msg = f"Known controls file {input.known_controls} not found. Skipping ROC plot for {wildcards.experiment}/{wildcards.contrast}."
            print(f"WARNING: {warn_msg}")
            with open(str(log), "w") as f:
                f.write(warn_msg)
            # touch(output.html_plot)
            sys.exit(0)

        print(
            f"Plotting ROC curve for {wildcards.experiment} - Contrast: {wildcards.contrast} (PLACEHOLDER)"
        )
        try:
            Path(output.html_plot).parent.mkdir(parents=True, exist_ok=True)
            print(
                "--> Placeholder rule: plot_roc_curve() function would be called here."
            )
            # Placeholder: Just log completion, no touch needed
            # touch(output.html_plot)
            print(f"ROC curve plot logged completion (placeholder): {output.html_plot}")

        except NameError:
            error_msg = "ERROR: plot_roc_curve function not found/imported (placeholder assumes it exists)."
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            # No touch needed for placeholder
            # touch(output.html_plot)
        except Exception as e:
            error_msg = f"Error plotting ROC curve for {wildcards.experiment}/{wildcards.contrast}: {e}"
            print(f"ERROR: {error_msg}")
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
            print(info_msg)
            with open(str(log), "w") as f:
                f.write(info_msg)
            # touch(output.report)
            sys.exit(0)

        print(f"Generating QC report for {wildcards.experiment} (PLACEHOLDER)")
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
            print(f"QC report logged completion (placeholder): {output.report}")

        except NameError:
            error_msg = "ERROR: generate_qc_report function not found/imported (placeholder assumes it exists)."
            print(error_msg)
            with open(str(log), "w") as f:
                f.write(error_msg)
            # No touch needed
            # touch(output.report)
        except Exception as e:
            error_msg = f"Error generating QC report for {wildcards.experiment}: {e}"
            print(f"ERROR: {error_msg}")
            with open(str(log), "w") as f:
                f.write(error_msg)
            # touch(output.report)
            raise e


# --- Rule All: Define Final Output Files ---


def get_final_outputs():
    """Dynamically collects all expected final output files based on config and targets."""
    final_files = []

    for experiment in TARGET_EXPERIMENTS:
        validation_info = get_validation_info(experiment)
        if validation_info["status"] == "failed":
            print(
                f"Skipping final output collection for failed experiment: {experiment}"
            )
            continue

        contrasts = get_contrast_names(SimpleNamespace(experiment=experiment))
        if (
            not isinstance(contrasts, list) or "error" in contrasts[0]
        ):  # Handle error case from get_contrast_names
            print(
                f"Skipping contrast-specific outputs for {experiment} due to contrast parsing issues."
            )
            contrasts = []

        # 1. Converted Analysis Results (per contrast)
        for contrast in contrasts:
            if not config.get("skip_rra", False):
                final_files.append(
                    OUTPUT_DIR / experiment / contrast / f"{contrast}_gMGK.csv"
                )

            # MLE requires design matrix AND skip flag check
            design_matrix_path = OUTPUT_DIR / experiment / "design_matrix.txt"
            if (
                not config.get("skip_mle", False)
                and design_matrix_path.exists()
                and design_matrix_path.stat().st_size > 0
            ):
                final_files.append(
                    OUTPUT_DIR / experiment / contrast / f"{contrast}_gMLE.csv"
                )

            if not config.get("skip_drugz", False):
                final_files.append(
                    OUTPUT_DIR / experiment / contrast / f"{contrast}_gDZ.csv"
                )

        # 2. QC Files (per experiment / per contrast)
        if not config.get("skip_qc", False):
            # Per-experiment QC plots
            final_files.append(
                OUTPUT_DIR / experiment / "qc" / f"{experiment}_sgrna_distribution.html"
            )
            final_files.append(
                OUTPUT_DIR / experiment / "qc" / f"{experiment}_gene_distribution.html"
            )
            final_files.append(
                OUTPUT_DIR / experiment / "qc" / f"{experiment}_gini_index.html"
            )
            final_files.append(
                OUTPUT_DIR / experiment / "qc" / f"{experiment}_qc_report.html"
            )  # Aggregate report

            # Per-contrast QC plots (ROC needs controls check)
            known_controls_path = BASE_DIR / experiment / "known_controls.csv"
            for contrast in contrasts:
                if known_controls_path.exists():  # Only add ROC if controls file exists
                    final_files.append(
                        OUTPUT_DIR
                        / experiment
                        / "qc"
                        / f"{experiment}_{contrast}_roc.html"
                    )

            # Per-sample FastQC reports (if FASTQ input)
            if validation_info.get("data_type") == "fastq":
                for sample in get_fastq_basenames(experiment):
                    final_files.append(
                        OUTPUT_DIR / experiment / "qc" / f"{sample}_fastqc.html"
                    )
                    # Add .zip as well if FastQC rule generates it
                    final_files.append(
                        OUTPUT_DIR / experiment / "qc" / f"{sample}_fastqc.zip"
                    )

    # Convert Path objects to strings for Snakemake input list
    return [str(f) for f in final_files]


# --- Final Target Rule ---
# Remove the temporary get_first_experiment_qc_targets function
# def get_first_experiment_qc_targets(): ...


rule all:
    input:
        get_final_outputs(),  # Use the function to gather all targets
    output:
        # Single flag file indicating completion of requested targets
        OUTPUT_DIR / "pipeline_complete.flag",
    run:
        print(f"CRISPR Analysis Pipeline Workflow Completed for targets:")
        for f in input:
            print(f"- {f}")
        print(f"Completion flag: {output}")
