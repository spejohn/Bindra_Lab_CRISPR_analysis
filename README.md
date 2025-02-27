# CRISPR Analysis Pipeline

A modular pipeline for CRISPR screening data analysis, supporting MAGeCK (RRA and MLE methods) and DrugZ analysis.

## Features

- Individual sample processing with Docker for better parallelization
- Comprehensive quality control and visualization
- Robust error handling and diagnostics
- Support for both read count tables and raw FASTQ files
- Automatic handling of dual input types (FASTQ and count files) with separate output directories
- Flexible input formats with automatic validation and fixing
- Conditional MAGeCK MLE analysis based on design matrix availability
- Standardized output directory structure with flattened organization
- Automatic CSV conversion of all analysis results for easier downstream processing
- **NEW: Snakemake workflow management for reproducible analyses**
- **NEW: Cross-platform compatibility (Windows, macOS, Linux)**
- **NEW: Enhanced Docker volume handling and path resolution**

## Directory Structure

The pipeline is organized into the following modules:

- `core/` - Core pipeline functionality
  - `config.py` - Configuration and constants
  - `logging_setup.py` - Logging configuration
  - `validation.py` - Input validation functions
  - `file_handling.py` - File and directory management utilities

- `docker/` - Docker-related functionality
  - `docker_utils.py` - Docker container management utilities
  - `mageck_docker.py` - MAGeCK Docker functions
  - `drugz_utils.py` - DrugZ execution utilities

- `analysis/` - Analysis functions
  - `sample_processing.py` - Sample processing and count aggregation
  - `mageck_analysis.py` - MAGeCK RRA and MLE analysis

- `qc/` - Quality control components
  - `plot_quality.py` - QA/QC plotting functions
  - `correlation_plots.py` - Guide- and gene-level correlation plots
  - `statistics.py` - Statistical analysis utilities

- `utils/` - General utility functions
  - `sample_utils.py` - Sample sheet generation and management
  - `diagnostic_utils.py` - Diagnostic and troubleshooting tools

## Output Directory Structure

The pipeline creates a standardized directory structure with flattened organization for better clarity:

```
<output_dir>/
├── <experiment_name>/                      # Main experiment directory (or <experiment_name>_FASTQ/ or <experiment_name>_RC/)
│   ├── <experiment_name>.count             # Merged count file for the entire experiment
│   ├── <experiment_name>.log               # Log file for the entire experiment
│   ├── <experiment_name>_samples.txt       # Sample sheet for the experiment
│   │
│   ├── <contrast_name1>/                   # First contrast directory
│   │   ├── <contrast_name1>_RRA.gene_summary.txt     # MAGeCK RRA results
│   │   ├── <contrast_name1>_gMGK.csv                 # MAGeCK RRA results in CSV format
│   │   ├── <contrast_name1>_MLE.gene_summary.txt     # MAGeCK MLE results (when applicable)
│   │   ├── <contrast_name1>_gMLE.csv                 # MAGeCK MLE results in CSV format
│   │   ├── <contrast_name1>_DrugZ.txt                # DrugZ results
│   │   └── <contrast_name1>_gDZ.csv                  # DrugZ results in CSV format
│   │
│   ├── <contrast_name2>/                   # Second contrast directory
│   │   ├── <contrast_name2>_RRA.gene_summary.txt
│   │   ├── <contrast_name2>_gMGK.csv
│   │   └── ...
│   │
│   └── qc/                                 # Quality control plots and reports
```

When processing both FASTQ and count files, separate experiment directories are created:

```
<output_dir>/
├── <experiment_name>_FASTQ/                # Results from FASTQ file analysis
│   ├── <experiment_name>_FASTQ.count
│   ├── <experiment_name>_FASTQ.log
│   ├── <experiment_name>_FASTQ_samples.txt
│   └── ...
│
├── <experiment_name>_RC/                   # Results from read count file analysis
│   ├── <experiment_name>_RC.count
│   ├── <experiment_name>_RC.log
│   ├── <experiment_name>_RC_samples.txt
│   └── ...
```

### File Naming Conventions

- **Log files**: `<experiment_name>.log`
- **Sample sheets**: `<experiment_name>_samples.txt`
- **Count files**: `<experiment_name>.count`
- **RRA results**: `<contrast_name>_RRA.gene_summary.txt` and `<contrast_name>_gMGK.csv`
- **MLE results**: `<contrast_name>_MLE.gene_summary.txt` and `<contrast_name>_gMLE.csv`
- **DrugZ results**: `<contrast_name>_DrugZ.txt` and `<contrast_name>_gDZ.csv`

## Input Directory Structure

The pipeline now supports a more organized input structure with experiment-specific subdirectories:

```
<input_dir>/
├── <experiment_name1>/           # First experiment directory
│   ├── library.csv               # CRISPR library file for this experiment
│   ├── contrasts.csv             # Contrast definitions file for this experiment
│   ├── design_matrix.txt         # Design matrix for MLE analysis (optional)
│   │
│   ├── fastq/                    # Directory containing FASTQ files
│   │   ├── sample1_R1.fastq.gz   # FASTQ files should follow standard naming
│   │   ├── sample2_R1.fastq.gz
│   │   └── ...
│   │
│   └── counts/                   # Directory containing count files
│       ├── other_sample1.count   # Count files for this experiment
│       ├── other_sample2.count
│       └── ...
│
├── <experiment_name2>/           # Second experiment directory
│   ├── library.csv               # CRISPR library file for this experiment
│   ├── contrasts.csv             # Contrast definitions file for this experiment
│   ├── design_matrix.txt         # Design matrix for MLE analysis (optional)
│   │
│   ├── fastq/                    # Directory containing FASTQ files
│   │   ├── sample1_R1.fastq.gz
│   │   └── ...
│   │
│   └── counts/                   # Directory containing count files
│       ├── other_sample1.count
│       └── ...
```

**Note**: The pipeline only recognizes "fastq" and "counts" as valid data directories. Generic directory names like "data" are not supported to avoid ambiguity.

When running the pipeline, the `experiment_name` parameter specifies which subdirectory to use for analysis:

```bash
python pipeline.py --input-dir /path/to/input --output-dir /path/to/output --experiment-name experiment1
```

This will:
1. Look for files in `/path/to/input/experiment1/`
2. Create output in `/path/to/output/experiment1/`

### Legacy Input Structures

The pipeline still supports the following legacy input structures for backward compatibility:

#### Option 1: Raw FASTQ files

```
<input_dir>/
├── library.csv                  # CRISPR library file (required)
├── contrasts.csv                # Contrast definitions file (required)
├── design_matrix.txt            # Design matrix for MLE analysis (optional)
├── fastq/                       # Directory containing FASTQ files
│   ├── sample1_R1.fastq.gz      # FASTQ files should follow standard naming
│   ├── sample2_R1.fastq.gz
│   └── ...
```

#### Option 2: Pre-processed count files

```
<input_dir>/
├── library.csv                  # CRISPR library file (required)
├── contrasts.csv                # Contrast definitions file (required)
├── design_matrix.txt            # Design matrix for MLE analysis (optional)
├── counts/                      # Directory containing count files
│   ├── sample1.count            # Count files should be tab-delimited
│   ├── sample2.count
│   └── ...
```

#### Option 3: Mixed input (both FASTQ and count files)

```
<input_dir>/
├── library.csv                  # CRISPR library file (required)
├── contrasts.csv                # Contrast definitions file (required)
├── design_matrix.txt            # Design matrix for MLE analysis (optional)
├── fastq/                       # Directory containing FASTQ files
│   ├── sample1_R1.fastq.gz
│   └── ...
├── counts/                      # Directory containing count files
│   ├── other_sample1.count
│   └── ...
```

### Required Input Files

1. **Library File (library.csv)**:
   - CSV format with columns for guide ID, gene, and sequence
   - Example:
     ```
     sgRNA,Gene,Sequence
     guide1,GENE1,ACGTACGTACGTACGTACGT
     guide2,GENE2,GTACGTACGTACGTACGTAC
     ```

2. **Contrasts File (contrasts.csv)**:
   - CSV format defining experimental contrasts
   - Must contain columns for contrast name, control samples, and treatment samples
   - Example:
     ```
     contrast,control,treatment
     drug_vs_control,control1|control2,drug1|drug2
     knockout_vs_wildtype,wt1|wt2,ko1|ko2
     ```

3. **Design Matrix (design_matrix.txt, optional)**:
   - Tab-delimited file for MLE analysis
   - First column must be 'Sample' followed by condition columns
   - Example:
     ```
     Sample  Treatment  Replicate
     control1    0    1
     control2    0    2
     drug1    1    1
     drug2    1    2
     ```

## Main Entrypoints

- `pipeline.py` - Main pipeline script
- `run_individual_sample.py` - Process a single sample (for parallelization)
- `run_qc.py` - Standalone QC analysis

## Usage

### Snakemake Runner

The simplest way to run the pipeline is using the Snakemake runner:

```bash
python analysis_pipeline/run_snakemake.py /path/to/screens_directory
```

Additional options:
- `-o, --output-dir`: Set custom output directory (default: same level as input directory, named "crispr_analysis_pipeline_results")
- `-j, --cores`: Number of CPU cores to use (default: auto-detected maximum cores)
- `--dryrun`: Show what would be done without executing
- `--skip-drugz`: Skip DrugZ analysis
- `--skip-qc`: Skip quality control checks
- `--skip-mle`: Skip MAGeCK MLE analysis
- `--configfile`: Use a custom configuration file

#### Automatic CPU Core Detection

The pipeline automatically detects and uses all available CPU cores on your system. You can:

1. Let the system automatically detect and use all available cores (default)
2. Specify a specific number with the `-j` option (e.g., `-j 4`)
3. Set the `CRISPR_MAX_CORES` environment variable to control maximum core usage:

```bash
# Set environment variable to use 8 cores maximum
export CRISPR_MAX_CORES=8

# Run pipeline (will use up to 8 cores)
python analysis_pipeline/run_snakemake.py /path/to/screens_directory
```

This is particularly useful for HPC environments where you may want to use all allocated cores for your job.

The pipeline automatically handles workflow locks, ensuring you can rerun the analysis anytime without worrying about lockfile errors that commonly occur if a previous run was interrupted.

```bash
# Run full pipeline
python pipeline.py --input-dir /path/to/input --output-dir /path/to/output --library-file /path/to/library.csv --experiment-name my_experiment

# Run pipeline with specific contrasts file
python pipeline.py --input-dir /path/to/input --output-dir /path/to/output --library-file /path/to/library.csv --contrasts-file /path/to/contrasts.csv --experiment-name my_experiment

# Run pipeline without DrugZ analysis
python pipeline.py --input-dir /path/to/input --output-dir /path/to/output --library-file /path/to/library.csv --experiment-name my_experiment --skip-drugz

# Run pipeline without MLE analysis
python pipeline.py --input-dir /path/to/input --output-dir /path/to/output --library-file /path/to/library.csv --experiment-name my_experiment --skip-mle

# Process individual sample
python run_individual_sample.py --fastq /path/to/sample.fastq --library-file /path/to/library.csv --output-dir /path/to/output

# Run QC on existing results
python run_qc.py --input /path/to/results --output /path/to/qc_output
```

### Example Workflow for Dual Input Types

If you have both FASTQ files and pre-processed count files in your input directory:

```bash
python pipeline.py --input-dir /path/to/mixed_input --output-dir /path/to/results --library-file /path/to/library.csv --experiment-name my_experiment
```

This will:
1. Detect both FASTQ files and count files in the input directory
2. Create two separate experiment directories:
   - `/path/to/results/my_experiment_FASTQ/` - For FASTQ-based analysis
   - `/path/to/results/my_experiment_RC/` - For read count-based analysis
3. Run the full analysis pipeline on both data types independently
4. Maintain separate log files and results for each analysis

### Example Workflow

1. **Prepare your input files**:
   - Place your FASTQ files in the input directory or ensure count files are organized
   - Create a library.csv file with guide information
   - Create a contrasts.csv file to define your experimental contrasts
   - (Optional) Create a design_matrix.txt for MLE analysis

2. **Run the pipeline**:
   ```bash
   python pipeline.py --input-dir /path/to/fastq_dir --output-dir /path/to/results --library-file /path/to/library.csv --experiment-name my_experiment
   ```

3. **Analyze the results**:
   - View the merged count table at `<output_dir>/<experiment_name>/<experiment_name>.count`
   - Check analysis results for each contrast in their respective directories
   - Explore the CSV files (`*_gMGK.csv`, `*_gMLE.csv`, `*_gDZ.csv`) for easy import into other tools
   - Review the log file at `<output_dir>/<experiment_name>/<experiment_name>.log`

## Understanding Analysis Results

The pipeline produces several types of analysis files that are saved in standardized formats:

### MAGeCK RRA Results (`*_RRA.gene_summary.txt` and `*_gMGK.csv`)

MAGeCK RRA (Robust Rank Aggregation) performs a non-parametric analysis to identify genes with significant guide enrichment or depletion:

- **Key columns**:
  - `Gene`: Gene identifier
  - `neg|score`: Depletion score (lower is more depleted)
  - `neg|lfc`: Depletion log fold change
  - `neg|p-value` and `neg|fdr`: Statistical significance for depletion
  - `pos|score`: Enrichment score (lower is more enriched)
  - `pos|lfc`: Enrichment log fold change
  - `pos|p-value` and `pos|fdr`: Statistical significance for enrichment

- **CSV format**:
  The `*_gMGK.csv` file contains the same information reformatted for easy import into spreadsheet software or R.

### MAGeCK MLE Results (`*_MLE.gene_summary.txt` and `*_gMLE.csv`)

MAGeCK MLE (Maximum Likelihood Estimation) uses a negative binomial model to estimate gene effects:

- **Key columns**:
  - `Gene`: Gene identifier
  - `Treatment|beta`: Effect size for each treatment/condition
  - `Treatment|wald-p-value` and `Treatment|fdr`: Statistical significance for each condition
  - `LRT p-value` and `LRT FDR`: Overall likelihood ratio test results

- **CSV format**:
  The `*_gMLE.csv` file contains the same information reformatted for easier import and analysis.

### DrugZ Results (`*_DrugZ.txt` and `*_gDZ.csv`)

DrugZ is specialized for analyzing drug-gene interactions in CRISPR screens:

- **Key columns**:
  - `GENE`: Gene identifier
  - `NORM_CD`: Normalized change in abundance
  - `Z_SCORE`: Z-score of the effect
  - `P_VALUE` and `FDR`: Statistical significance

- **CSV format**:
  The `*_gDZ.csv` file contains the same information with standardized column headers for consistency.

### Count Table (`<experiment_name>.count`)

The merged count table is a tab-delimited file containing read counts for each guide across all samples:

- **Format**:
  - First column: Guide IDs
  - Second column: Gene symbols
  - Remaining columns: Read counts for each sample

### Sample Sheet (`<experiment_name>_samples.txt`)

The sample sheet lists all samples included in the analysis with their file paths:

- **Format**:
  - `sample`: Sample name
  - `fastq`: Path to the FASTQ file or count file
  - `condition`: Experimental condition (if available)

### Log File (`<experiment_name>.log`)

The log file contains a detailed record of the analysis process, including:

- Command-line parameters used
- Processing steps and timestamps
- Warnings and errors encountered
- Summary statistics for each analysis

## MAGeCK MLE Analysis

The pipeline automatically detects design matrix files in the same directory as contrast tables. If a design matrix is found, it will run MAGeCK MLE analysis in addition to the standard RRA analysis.

Requirements for design matrix files:
- Must be named following patterns like `design_matrix.txt`, `design.txt`, etc.
- Must be in the same directory as the contrast table
- Must follow the MAGeCK MLE design matrix format with columns for Sample and condition variables

## Multi-Screen Analysis

The pipeline now supports automatic detection and analysis of multiple screens:

```bash
# Run analysis on all screens in a directory
python analysis_pipeline/run_multiple_screens.py /path/to/screens_directory
```

This will:
1. Automatically detect all directories containing CRISPR screen data
2. Find library and contrast files for each screen
3. Determine which screens have already been analyzed
4. Only process screens that need analysis
5. Generate a summary report of all screens

## Smart Path Detection

The pipeline can now automatically find necessary files:

```bash
# Simplified usage with auto-detection of outputs and library
python analysis_pipeline/pipeline.py /path/to/screen_directory
```

The pipeline will:
- Use `<input_dir>/library.csv` as the default library file if not specified
- Use `<input_dir>/contrasts.csv` as the default contrasts file if not specified
- Create output in `<input_dir>/results` by default
- Detect if the screen has already been analyzed and skip if appropriate

## Cross-Platform Compatibility

This pipeline is designed to work seamlessly across Windows, macOS, and Linux operating systems. Special attention has been given to path handling to ensure compatibility with Docker volumes on all platforms.

### Windows Path Handling for Docker

When running on Windows, the pipeline automatically converts Windows-style paths to Docker-compatible paths:

- Windows path: `C:\Users\username\data`
- Docker path: `/c/Users/username/data`

This conversion happens transparently in the background, so you can use regular Windows paths throughout the pipeline.

### Path Conversion Function

The pipeline provides a utility function for path conversion in `core/utils.py`:

```python
from analysis_pipeline.core.utils import convert_win_path_to_docker

# Convert a Windows path to Docker-compatible format
windows_path = r"C:\Users\username\data"
docker_path = convert_win_path_to_docker(windows_path)
# Result: "/c/Users/username/data"
```

### Testing Path Conversion

You can verify path conversion is working correctly by running the included test script:

```bash
python analysis_pipeline/test_docker_paths.py
```

This will output the conversion results for various test paths and your current working directory.

### Docker Volume Mounting

The pipeline handles Docker volume mounting differently depending on the operating system:

- **Windows**: Automatically converts paths for Docker compatibility
- **macOS/Linux**: Uses paths directly with minimal modification

This ensures that file access within Docker containers works correctly regardless of the host operating system.

## Requirements

- Python 3.8+
- Docker
- Required Python packages are listed in requirements.txt

## Docker Containers

The pipeline uses containerized execution of analysis tools, which provides several advantages:

- Consistent execution environment across different systems
- No need to install tools and dependencies locally
- Easier deployment and reproducibility

### Available Containers

- **MAGeCK Container**: For MAGeCK RRA and MLE analyses
- **DrugZ Container**: For DrugZ analysis

### Building Docker Images

To build the Docker images, run:

```bash
python analysis_pipeline/docker/build_images.py
```

Or manually:

```bash
cd analysis_pipeline/docker
docker-compose build
```

### Docker Usage

Docker is required for the pipeline to function properly as all analysis tools run in Docker containers. The pipeline automatically uses Docker for all analyses.

```bash
python analysis_pipeline/pipeline.py \
    --input-dir /path/to/fastq \
    --output-dir /path/to/output \
    --library-file /path/to/library.csv \
    --experiment-name my_experiment
```

You can also use the DrugZ Docker container directly:

```bash
python analysis_pipeline/analysis/run_drugz_docker.py \
    --input /path/to/count_file.csv \
    --design /path/to/design_file.txt \
    --output /path/to/output \
    --prefix my_experiment
```

### Docker Directory Structure

```
analysis_pipeline/docker/
├── docker-compose.yml           # Docker Compose configuration
├── build_images.py              # Script to build Docker images
└── dockerfiles/
    ├── Dockerfile.mageck        # Dockerfile for MAGeCK
    └── Dockerfile.drugz         # Dockerfile for DrugZ
```

## Snakemake Workflow

The pipeline includes a Snakemake workflow for improved reproducibility and dependency management:

```bash
# Run the workflow with default settings
python analysis_pipeline/run_snakemake.py /path/to/screens_directory

# Run with multiple cores
python analysis_pipeline/run_snakemake.py /path/to/screens_directory -j 8

# Do a dry run to see what would be executed
python analysis_pipeline/run_snakemake.py /path/to/screens_directory --dryrun

# Specify custom output directory
python analysis_pipeline/run_snakemake.py /path/to/screens_directory -o /path/to/output
```

### Advantages of Snakemake

- **Automatic dependency tracking**: Only runs necessary steps based on file timestamps
- **Parallel execution**: Process multiple samples simultaneously
- **Resume functionality**: Pick up where you left off if a run is interrupted
- **Workflow visualization**: Generate DAG visualizations of your workflow
- **Reproducibility**: Ensures the same steps are run in the same order

### Configuration

The workflow can be configured using a `config.yaml` file or command-line parameters:

```yaml
# Example config.yaml
input_dir: "input"
output_dir: "crispr_analysis_pipeline_results"
threads: 8
skip_drugz: false
# Docker is always enabled as it's required for analysis
```

To use a custom config file:

```bash
python analysis_pipeline/run_snakemake.py /path/to/screens_directory --configfile /path/to/config.yaml
```

## Installation

### Dependencies

Install the required dependencies:

```bash
pip install -r requirements.txt
```

For Snakemake to work properly on Windows, you may need to install dependencies using Conda:

```bash
conda install -c bioconda -c conda-forge snakemake
```

### Windows Compatibility

The pipeline now includes enhanced Windows support:

- **Path Handling**: Automatic conversion of Windows paths to Docker-compatible format
- **Volume Mounting**: Improved Docker volume mounting for Windows host paths
- **Debugging**: Additional logging for Docker command execution and path conversion
- **Docker Desktop**: Requires Docker Desktop on Windows with WSL 2 backend for optimal performance

To verify Docker is properly configured on Windows:

```bash
python -c "from analysis_pipeline.docker.docker_utils import verify_docker; verify_docker()"
```

### Note on Linter Errors in Snakefile

The Snakefile may show linter errors in some IDEs due to Snakemake's unique syntax, which extends Python. These errors can be safely ignored as they don't affect the functionality of the workflow. Common false positive errors include:

- Triple-quoted shell commands ("""...""")
- The use of rule blocks and wildcards
- Snakemake-specific directives like configfile, rule, and input/output sections

If modifying the Snakefile, it's recommended to test changes with a dry run:

```bash
python analysis_pipeline/run_snakemake.py --dryrun /path/to/screen_directory
```

### Workflow Recovery and Locks

The pipeline will automatically handle workflow locks that might occur if a previous run was interrupted unexpectedly (due to power loss, system crash, etc.). This eliminates the need for manual intervention to unlock workflows.

If you encounter any unusual behavior after a workflow was interrupted, you can still perform a clean restart by:

1. Running with `--dryrun` to see the planned workflow
2. If needed, manually removing the `.snakemake` directory in your output folder
   
This automatic lock handling makes the pipeline more robust and user-friendly, especially for those new to Snakemake workflows.