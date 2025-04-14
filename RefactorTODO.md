# CRISPR Analysis Pipeline Refactoring Plan

This document outlines the plan to refactor the CRISPR analysis pipeline codebase to support two distinct execution modes: a Snakemake-based mode for HPC/batch processing and a standalone Python mode for local/targeted analysis.

## Target Directory Structure

This illustrates the expected organization of input data:

```
<Main_dir>/                     # Base directory (Input for run_snakemake.py)
├── <First_Last_LibraryID_1>/   # Specific experiment directory (Input for run_crispr_analysis.py)
│   ├── library.csv             # Required if fastq/ exists (for gene mapping)
│   ├── contrasts.csv           # Required (defines comparisons)
│   ├── design_matrix.txt       # Optional (required for MAGeCK MLE)
│   ├── fastq/                  # Optional directory for FASTQ files
│   │   ├── sample1_R1.fastq.gz
│   │   ├── sample1_R2.fastq.gz # Optional (paired-end)
│   │   └── ...
│   └── experiment_counts.csv   # Optional (alternative to fastq/; single file)
│
├── <First_Last_LibraryID_2>/   # Another experiment directory...
│   ├── library.csv
│   ├── contrasts.csv
│   └── experiment_counts.csv
│
└── ...
```

**Notes:**
- `run_snakemake.py` takes `<Main_dir>` and processes subdirectories.
- `run_crispr_analysis.py` takes one or more `<First_Last_LibraryID>` paths.
- Each experiment directory can contain *either* a `fastq/` directory *or* an `experiment_counts.csv` file, or potentially both if different analyses are intended (though the exact handling of both needs definition).
- File names (`library.csv`, `contrasts.csv`, `design_matrix.txt`, `experiment_counts.csv`) should be consistent or configurable.

## Goal Architecture

1.  **Core Analysis Module(s):** Python functions encapsulating individual analysis steps (FastQC, MAGeCK count, MAGeCK RRA/MLE, DrugZ, etc.), callable independently from Snakemake, potentially using `subprocess` and Docker.
2.  **`Snakefile`:** Defines the end-to-end workflow using the core analysis functions. Handles experiment discovery (within a base directory), optional filtering based on a target list, dependencies, and parallelism.
3.  **`run_snakemake.py` (HPC/Batch Entry Point):**
    *   Requires Snakemake installation.
    *   Input: `Main_dir`, optional list of specific screen names (`--screens`), output dir, cores, skips, overwrite flags.
    *   Executes `snakemake` command, passing configuration to the `Snakefile`.
4.  **`run_crispr_analysis.py` (Local/Point-and-Shoot Entry Point):** (Refactor/rename `run_multiple_screens.py`)
    *   Does NOT require Snakemake installation.
    *   Input: One or more specific `<First_Last_LibraryID>` directory paths, output dir, overwrite flag, optional simple parallelism control.
    *   Directly calls the core analysis functions for each specified directory.

## Refactoring Plan Outline

### Phase 1: Develop Core Analysis Functions
- **Purpose:** Create modular, reusable Python functions for each distinct analysis step, callable by both Snakemake (Phase 2) and the standalone script (Phase 4).
- **Location:** Primarily within `core/analysis_steps.py` (or similar module).
- **General Design:** Each function should:
    - Accept specific input file paths and parameters.
    - Accept output file path(s) or directory.
    - Encapsulate the logic to run a specific tool (e.g., construct command line).
    - Utilize a helper function (`execute_in_container`) for running commands via Apptainer/Docker.
    - Perform basic validation of inputs/outputs.
    - Return status and/or path(s) to generated output(s).
- **Tasks:**
    - **Container Execution Helper:**
        - `execute_in_container(command_list, container_image, mount_dirs, working_dir)`: Helper function to run a command list within a specified Apptainer/Docker container, handling necessary volume mounts and working directory. Uses `subprocess.run`. Should handle logging of command execution and errors. Indirectly used by Phase 2 & 4.
    - **Input Validation & Parsing:**
        - `validate_experiment_structure(experiment_dir)`: Renamed function. Checks for the *existence* and 
        *basic structure* of required files based on context (`contrasts.csv` always; `library.csv` if 
        `fastq/` present; `*_rc.csv` if `fastq/` absent; `design_matrix.csv` optional). Returns a dictionary 
        detailing found files, determined data type (fastq/rc), and potential for MLE, or raises validation 
        errors on structure (e.g., missing required files, or **if both fastq/ and *_rc.csv are present**). *Detailed content validation (e.g., column names/formats) will occur in 
        subsequent functions after potential file conversions.* Leverages basic file system checks. Called by Phase 4 (local script) for upfront validation. Potentially called by Phase 2 (Snakemake) in an initial validation rule.
        - `parse_contrasts(contrasts_path)`: Reads the *converted* tab-delimited contrasts file. Returns a 
        list of tuples/dictionaries representing each contrast. Called by Phase 2 (Snakemake rules for diff. analysis) & Phase 4 (local script) before running diff. analysis.
    - **FastQ Processing Functions:** (Called by Phase 2 rules & Phase 4 script when processing FASTQ)
        - `run_fastqc(fastq_path, output_dir)`: Takes path to a single `.fastq.gz` file. Constructs and executes FastQC command via `execute_in_container`. Requires appropriate FastQC container image.
        - `run_mageck_count(r1_fastq, library_path, output_prefix, r2_fastq=None, sample_name=None)`: Takes paths for R1 (required), library (required), and optional R2. Constructs and executes `mageck count` command via `execute_in_container`. Requires MAGeCK container image. Handles gzipped input.
    - **Count Data Processing Functions:** (Called by Phase 2 rule & Phase 4 script after sample counts generated)
        - `aggregate_mageck_counts(sample_count_files, output_path)`: Takes a list of paths to individual MAGeCK count outputs. Merges them into a single table compatible with downstream tools. Writes to `output_path`.
    - **Differential Analysis Functions:** (Called by Phase 2 rules & Phase 4 script after counts available)
        - `run_mageck_rra(count_path, contrast_info, output_prefix)`: Takes path to aggregated counts, a single contrast dictionary (from `parse_contrasts`), and output prefix. Constructs and executes `mageck test` (RRA) command via `execute_in_container`.
        - `run_mageck_mle(count_path, design_matrix_path, contrast_info, output_prefix)`: Takes path to counts, design matrix, a contrast dictionary, and output prefix. Constructs and executes `mageck mle` command via `execute_in_container`.
        - `run_drugz(count_path, contrast_info, output_prefix)`: Takes path to counts, a contrast dictionary, and output prefix. Constructs and executes DrugZ command via `execute_in_container`. Requires DrugZ container image.
    - **Utility Functions:**
        - `convert_csv_to_tab(input_csv_path, output_txt_path)`: **[EXISTING]** Implemented in `core/file_handling.py::convert_file_to_tab_delimited`. Handles conversion of `contrasts.csv` and `design_matrix.csv` to tab-delimited `.txt` files required by MAGeCK/DrugZ. Called by Phase 2 (Snakemake rule) & Phase 4 (local script) early in the process for relevant files.
            - Implements whitespace cleaning on values.
            - Supports standard and multi-column contrast CSV formats (merging columns as needed).
        - `make_count_table(rc_file, output_dir)`: **[EXISTING]** Implemented in `core/file_handling.py::make_count_table`. Converts an input CSV count table (`*_rc.csv` or `experiment_counts.csv`) to the required tab-delimited `.count.txt` format. Called by Phase 2 rule.
    - **Result Conversion Functions:** (Called by Phase 2 conversion rules)
        - `convert_mageck_summary_to_csv(input_txt_path, output_csv_path, analysis_type)`: Parses MAGeCK RRA/MLE `.gene_summary.txt` and saves key columns to a standardized CSV format (`_gMGK.csv` or `_gMLE.csv`).
        - `convert_drugz_summary_to_csv(input_txt_path, output_csv_path)`: Parses DrugZ `.drugz.txt` output and saves key columns to a standardized CSV format (`_gDZ.csv`).
    - **QC & Reporting Functions:** (Called by Phase 2 rules & Phase 4 script, potentially skipped via flag)
        - `plot_sgRNA_distribution(count_path, output_html_path)`: Renamed. Generates interactive histogram/density plot (e.g., using Plotly) of read counts per guide from the aggregated count file. Saves as HTML.
        - `plot_gene_distribution(count_path, output_html_path)`: New. Generates interactive plot of read counts summed per gene from the aggregated count file. Saves as HTML.
        - `plot_gini_index(count_path, output_html_path)`: Calculates Gini index from counts and generates an interactive visualization (e.g., Lorenz curve). Saves as HTML.
        - `plot_roc_curve(mageck_results_path, known_controls_path, output_html_path)`: **[TODO]** Generates interactive ROC curve based on MAGeCK results and known controls (if provided). Saves as HTML. Requires scikit-learn.
        - `generate_qc_report(...)`: **[TODO]** Potentially a function to aggregate various QC metrics and plots into a single HTML report (e.g., using `plotly.io.to_html` fragments or a reporting library).

### Phase 2: Refactor `Snakefile`
- **Purpose:** Define the complete, reproducible workflow logic using Snakemake, calling core analysis functions where possible.
- **Tasks:**
    - **Experiment Discovery:** Use `glob_wildcards` to identify potential `experiment_dir` subdirectories within the configured `base_dir`.
    - **Target Filtering:** Implement logic to process only experiments listed in the `target_screens` config list, if provided. Otherwise, process all discovered experiments.
    - **Input Validation Rule/Function:** Create logic (e.g., an input function or initial checkpoint rule) that runs *per discovered experiment* to:
        - Verify `contrasts.csv` exists. Log error and skip experiment if not.
        - Check for presence of `fastq/` vs `experiment_counts.csv`.
        - If `fastq/` exists, verify `library.csv` exists. Log error and skip if not.
        - If `fastq/` does *not* exist, verify `experiment_counts.csv` exists. Log error and skip if not.
        - Check for `design_matrix.txt` and set flag/config for conditional MLE.
        - Ensure downstream rules depend on successful validation.
    - **FastQC Rule:** Define rule to run FastQC on each `.fastq.gz` file within the `fastq/` subdirectory. Tool must handle gzipped input directly or include an unzipping step. Output FastQC reports to the experiment's output directory.
    - **MAGeCK Count Rule:** Define rule per sample. Input: R1 `.fastq.gz` (and optionally R2), `library.csv`. Output: Sample-specific count file (e.g., `.count.txt`). Handle gzipped input.
    - **Aggregate Counts Rule:** Define rule per experiment. Input: All sample-specific count files for the experiment. Output: Merged count file (e.g., `<experiment_name>.count`).
    - **Analysis Rules (RRA, MLE, DrugZ):** Define rules per contrast found in `contrasts.csv`.
        - Input: Merged count file (or `experiment_counts.csv`), processed `contrasts.csv`.
        - MLE rule additionally requires `design_matrix.txt` and should be skipped if file is absent or `skip_mle` config is true.
        - RRA/DrugZ rules should be skipped if corresponding skip flags are true.
        - Output: Result files named appropriately (e.g., `<contrast>_gMGK.csv`, `<contrast>_gMLE.csv`, `<contrast>_gDZ.csv`).
    - **Result Conversion Rules:** Define rules that depend on the native text outputs of RRA, MLE, and DrugZ rules.
        - Calls the corresponding Phase 1 conversion functions (`convert_mageck_summary_to_csv`, `convert_drugz_summary_to_csv`).
        - Output: Standardized CSV files (`_gMGK.csv`, `_gMLE.csv`, `_gDZ.csv`) in the contrast-specific directory.
    - **QC Rules:** Define rules to generate QC plots/metrics, conditionally executed based on the `skip_qc` config flag.
        - Input: Aggregated count file, MAGeCK results, potentially known controls file.
        - Calls corresponding Phase 1 QC functions (`plot_sgRNA_distribution`, `plot_gene_distribution`, `plot_gini_index`, `plot_roc_curve`, etc.).
        - Output: Interactive HTML plot files (e.g., `*_sgrna_dist.html`, `*_gene_dist.html`, `*_gini.html`, `*_roc.html`) or metric files (.txt, .json) in the experiment's output directory.
    - **Container Integration:** Ensure rules executing external tools specify execution via Apptainer (`--use-apptainer` will be passed from wrapper).
    - **Modularity:** Use core analysis functions (from Phase 1) via `script:` directive where possible to keep rule definitions clean.

### Phase 3: Refactor `run_snakemake.py`
- **Purpose:** Provide a user-friendly command-line interface to configure and launch the Snakemake workflow, primarily for HPC/batch usage.
- **Tasks:**
    - **System Checks:** Implement initial checks for Apptainer (or Docker) availability before proceeding.
    - **Argument Parsing:** Update `argparse` to accept:
        - `main_dir` (positional argument, base directory containing experiments).
        - `--screens` (optional, `nargs='*'`, list of specific experiment directory names within `main_dir` to process).
        - `--output-dir` (optional, path for results).
        - `--cores` (optional, number of cores for Snakemake).
        - Skip flags (`--skip-qc`, `--skip-mle`, `--skip-drugz`).
        - `--overwrite`

### Phase 4: Refactor `Snakefile` for Native Containerization
- **Purpose:** Modify the `Snakefile` to use Snakemake's native container support (`container:` and `--use-apptainer`/`--use-docker`) instead of passing a custom flag to Python functions.
- **Tasks:**
    - **Define Container URIs:** Define container image URIs (e.g., `FASTQC_IMAGE = "docker://..."`) at the top of the `Snakefile` or import them from a config module.
    - **Remove Custom Flag:** Remove the `config.setdefault("use_apptainer", ...)` logic. Container runtime selection will be controlled by command-line flags passed to `snakemake`.
    - **Refactor Rules (FastQC, MAGeCK count, RRA, MLE, DrugZ):**
        - Add a `container:` directive specifying the appropriate image URI to each rule executing an external tool.
        - Remove the `use_apptainer=...` entry from the `params:` section.
        - Replace the `run:` block (which calls Python functions like `run_fastqc`) with a `shell:` block.
        - **Command Construction:** Move the logic for building the command line (including handling options like paired-end reads, library formats, contrast samples) from the `core.analysis_steps.py` functions into the `shell:` directive strings, using Snakemake's formatting capabilities (`{input}`, `{output}`, `{params}`, `{threads}`, `{wildcards}`).
        - **Helper Functions (Optional):** If command construction becomes too complex within the `shell:` directive, define Python helper functions *within the Snakefile* that return formatted command strings or parts of commands, and call these from the `shell:` directive (e.g., `shell: helper_function(wildcards, input, output)`).
    - **Address Conditional Execution (`skip_*` flags):**
        - Remove the `if config["skip_qc"]:` (and similar) checks from *within* the `run:` blocks (as these blocks are being replaced by `shell:`).
        - Ensure that workflow skipping (based on `skip_qc`, `skip_mle`, etc.) is handled by conditionally requesting output files in `rule all` (via the `get_final_outputs` function) based on the configuration. The `run_snakemake.py` script will pass these config flags to Snakemake.
    - **Update Workflow Execution:** Modify `run_snakemake.py` and/or `run_workflow.sh` to pass `--use-apptainer` or potentially `--use-docker` (as the default if neither is specified and containers are needed) to the `snakemake` command, instead of relying on the removed `--config use_apptainer=...`.

### Phase 5: Refactor/Create `run_crispr_analysis.py` (Local)
- Rename `run_multiple_screens.py` or create new script.
- Adapt argument parsing: Accept list of specific experiment directory paths.
- Implement main loop iterating through provided directories.
- Inside loop:
    - Validate input directories.
    - Check for existing output and handle `--overwrite`.
    - Determine required steps based on available files and skip flags.
    - Call core analysis functions (Phase 1) sequentially.
    - Manage output paths.
- Implement logging and summary.
- (Optional) Add simple parallelism (e.g., `multiprocessing.Pool`).

### Phase 6: Update Documentation & Testing
- Update `README.md`, `docs/INPUT_OUTPUT.md`, etc., to reflect the new architecture and usage.
- Develop/update unit/integration tests for core functions and both execution modes. 