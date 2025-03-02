# CRISPR Analysis Pipeline

A modular pipeline for CRISPR screening data analysis, supporting MAGeCK (RRA and MLE methods) and DrugZ analysis.

**Current Version: 0.1.0**

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
- **NEW: Improved parallel processing for both sample counting and contrast analysis**
- **NEW: Code cleanup and consistent function naming**
- **NEW: File format conversion utility for design matrices and contrast tables**

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

The pipeline supports multiple input directory structures, depending on your data type and analysis needs.

### File Requirements by Input Type

#### When Starting with FASTQ Files:
- **Library file is required** - Used to map sgRNA sequences to genes
- Contrast table - Defines the experimental comparisons
- Design matrix (optional) - For MAGeCK MLE analysis

#### When Starting with Count Files:
- **Library file is NOT required** - Mapping has already been completed
- The pipeline now explicitly checks for FASTQ files and only requires a library file when they are present
- Contrast table - Defines the experimental comparisons
- Design matrix (optional) - For MAGeCK MLE analysis

### Directory Structure Options

#### Option 1: FASTQ input

```
<input_dir>/
├── library.csv                  # CRISPR library file (REQUIRED for FASTQ analysis)
├── contrasts.csv                # Contrast definitions file (required)
├── design_matrix.txt            # Design matrix for MLE analysis (optional)
├── fastq/                       # Directory containing FASTQ files
│   ├── sample1_R1.fastq.gz      # Forward reads
│   ├── sample1_R2.fastq.gz      # Reverse reads (for paired-end data)
│   └── ...
```

#### Option 2: Count file input

```
<input_dir>/
├── contrasts.csv                # Contrast definitions file (required)
├── design_matrix.txt            # Design matrix for MLE analysis (optional)
├── counts/                      # Directory containing count files
│   ├── sample1.count
│   ├── sample2.count
│   └── ...
```

#### Option 3: Mixed input (both FASTQ and count files)

```
<input_dir>/
├── library.csv                  # CRISPR library file (REQUIRED only for FASTQ processing)
├── contrasts.csv                # Contrast definitions file (required)
├── design_matrix.txt            # Design matrix for MLE analysis (optional)
├── fastq/                       # Directory containing FASTQ files
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz      # For paired-end data
│   └── ...
├── counts/                      # Directory containing count files
│   ├── other_sample1.count
│   └── ...
```

### Required Input Files

1. **Library File (library.csv)** - REQUIRED ONLY when processing FASTQ files:
   - CSV format with columns for guide ID, gene, and sequence
   - Used to map sgRNA sequences to their target genes
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
   - **Note**: The pipeline expects tab-delimited text files for analysis. Use the file conversion utility to convert CSV files to the required format.

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
   - **Note**: The pipeline expects tab-delimited text files for analysis. Use the file conversion utility to convert CSV files to the required format.

## Main Entrypoints

- `pipeline.py` - Main pipeline script
- `run_individual_sample.py` - Process a single sample (for parallelization)
- `run_qc.py` - Standalone QC analysis
- `convert_input_files.py` - Utility to convert CSV files to tab-delimited format

## Usage

### Command Types and Parameter Requirements

The pipeline offers two types of commands:

1. **High-Level "Run Everything" Commands** - Automatically detect as much as possible:
   - `run_snakemake.py` - Just point to a directory, everything is auto-detected
   - `pipeline.py` (with no specific selections) - Auto-detects experiment names, contrasts, etc.

2. **Targeted/Specific Commands** - Require explicit selection parameters:
   - `run_individual_sample.py` - Needs explicit sample selection
   - Analysis variant commands - For specialized processing

### Snakemake Runner (Recommended)

The simplest way to run the pipeline is using the Snakemake runner, which automatically detects everything:

```bash
# Most basic form - just point to your screens directory
python analysis_pipeline/run_snakemake.py /path/to/screens_directory

# Optionally specify output directory and cores
python analysis_pipeline/run_snakemake.py /path/to/screens_directory -o /path/to/output -j 8
```

The runner automatically:
- Detects experiment directories (no need to specify experiment names)
- Identifies input types (FASTQ files, count files, or both)
- Locates contrast files and library files
- Adjusts requirements based on input (e.g., library file is only required for FASTQ processing)
- Runs appropriate analysis steps for each experiment

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

### Standard Pipeline Usage

```bash
# Most basic form - processes everything in specified directory
# (folder name is used as experiment name, files auto-detected)
python pipeline.py --input-dir /path/to/experiment_dir

# Specify output directory (input directory name still used as experiment name)
python pipeline.py --input-dir /path/to/experiment_dir --output-dir /path/to/output

# Override auto-detection with explicit experiment name
python pipeline.py --input-dir /path/to/input --output-dir /path/to/output --experiment-name my_experiment

# Specify a non-default contrasts file
python pipeline.py --input-dir /path/to/input --contrasts-file /path/to/custom/contrasts.csv

# Skip certain analysis types
python pipeline.py --input-dir /path/to/input --skip-drugz
python pipeline.py --input-dir /path/to/input --skip-mle

# Custom library file location (only needed if not in standard location)
python pipeline.py --input-dir /path/to/input --library-file /path/to/custom/location/library.csv
```

### Targeted Command Usage

For specific operations, use more detailed parameters:

```bash
# Process individual sample (requires explicit file paths)
python run_individual_sample.py --fastq /path/to/sample.fastq --library-file /path/to/library.csv --output-dir /path/to/output

# Run QC on existing results (requires specific input/output)
python run_qc.py --input /path/to/results --output /path/to/qc_output

# Run multiple specific screens (requires directory path)
python run_multiple_screens.py /path/to/screens_directory
```

### Example Workflow for Dual Input Types

If you have both FASTQ files and pre-processed count files in your input directory:

```bash
# Auto-detection approach (recommended)
python pipeline.py --input-dir /path/to/mixed_input

# With explicit output directory
python pipeline.py --input-dir /path/to/mixed_input --output-dir /path/to/results

# With explicit experiment name (optional override)
python pipeline.py --input-dir /path/to/mixed_input --output-dir /path/to/results --experiment-name my_experiment
```

This will:
1. Detect both FASTQ files and count files in the input directory
2. Create two separate experiment directories:
   - `/path/to/results/my_experiment_FASTQ/` - For FASTQ-based analysis
   - `/path/to/results/my_experiment_RC/` - For read count-based analysis
3. Run the full analysis pipeline on both data types independently
4. Maintain separate log files and results for each analysis

### File Detection and Error Handling

The pipeline automatically looks for required files in standard locations:
- `library.csv`: Required only when processing FASTQ files
- `contrasts.csv`: Required for all analyses
- `design_matrix.txt`: Optional, used for MLE analysis if present

If required files are missing, the pipeline will:
- Display a clear error message explaining which file is missing
- Indicate whether the file is required based on the detected input types
- Exit gracefully without starting the analysis

### Example Workflow

1. **Prepare your input files**:
   - Place your FASTQ files in the input directory or ensure count files are organized
   - Create a library.csv file with guide information (required only for FASTQ processing)
   - Create a contrasts.csv file to define your experimental contrasts
   - (Optional) Create a design_matrix.txt for MLE analysis

2. **Run the pipeline** (simple auto-detect approach):
   ```bash
   python pipeline.py --input-dir /path/to/experiment_dir
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

## Input File Formats

### Contrast Table Format

Contrast tables define pairwise comparisons between treatment and control groups and must follow this format:

- **File type**: Tab-separated text file (`.txt`), not CSV
- **Required columns**: `contrast`, `control`, `treatment` (tab-separated)
- **Multiple sample format**: Comma-separated (e.g., `sample1,sample2`) within each column
- **Example**:
  ```
  contrast	control	treatment
  test_contrast	control1,control2	treatment1,treatment2
  ```

This format is used by both standard MAGeCK analysis and DrugZ analysis. The file must be tab-delimited with each column separated by a tab character, not a comma.

### Design Matrix Format

Design matrices are used specifically by MAGeCK MLE for more complex experimental designs:

- **File type**: Tab-separated text file (`.txt`)
- **Required format**: A column named `Sample` followed by condition columns
- **Each row**: Represents one sample and its associated condition values
- **No sample grouping**: Unlike contrast tables, each sample gets its own row
- **Example**:
  ```
  Sample	condition
  control1	control
  control2	control
  treatment1	treatment
  treatment2	treatment
  ```

The design matrix allows MAGeCK to model more complex experimental variables beyond simple treatment vs. control comparisons.

### File Naming Requirements

- **Contrast tables**: Should be named following patterns like `contrasts.txt`, `contrast_table.txt`
- **Design matrices**: Should be named following patterns like `design_matrix.txt`, `design.txt`

These files should be placed in the experiment's input directory for automatic detection by the pipeline.

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

### Input Requirements for Snakemake

The workflow automatically detects the type of input in each experiment directory and adjusts its requirements accordingly:

1. **For experiments with FASTQ files**:
   - **Library file is REQUIRED** - The `library.csv` file must be present in the experiment directory
   - The workflow will run the complete pipeline including the counting step

2. **For experiments with count files only**:
   - **Library file is NOT required** - The workflow will skip the counting step
   - The workflow will proceed directly to MAGeCK test and DrugZ analysis

3. **For experiments with both FASTQ and count files**:
   - **Library file is REQUIRED** for processing the FASTQ files
   - The workflow will run the appropriate steps for each input type

For detailed directory structure examples, see the [Input Directory Structure](#input-directory-structure) section.

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

## Recent Updates

### Code Cleanup and Function Naming

The codebase has been cleaned up to improve maintainability:

- Standardized function names and eliminated redundant code
- Added proper version tracking (starting with version 0.1.0)
- Improved consistency in function interfaces and parameter names

### Parallel Processing Improvements

- Sample counting now supports parallel processing with configurable worker count
- Contrast analysis uses parallel processing for faster execution
- Automatic detection of Snakemake environment to disable internal parallelism when appropriate

### Error Handling and Recovery

- Enhanced error recovery with retry mechanisms for file operations
- Improved logging with progress reporting
- Better handling of edge cases in file detection and processing

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Utilities

### File Format Conversion

The pipeline includes a utility script to convert CSV-formatted design matrices and contrast tables to the tab-delimited format required by MAGeCK and DrugZ:

```bash
# Convert a single file
python convert_input_files.py --file /path/to/your/design_matrix.csv

# Convert all files in a directory
python convert_input_files.py --dir /path/to/your/input_directory

# Specify an output directory
python convert_input_files.py --file /path/to/your/contrast_table.csv --output /path/to/output/directory
```

For more information, see the [File Conversion Utility Documentation](README_file_conversion.md).

## Working with Paired-End Reads

The pipeline supports paired-end sequencing data, which is common in many CRISPR screening experiments:

### How Paired-End Reads Are Processed

1. **FASTQ File Naming**:
   - The pipeline identifies paired files based on naming conventions
   - Forward reads typically end with `_R1.fastq.gz` (or similar)
   - Reverse reads typically end with `_R2.fastq.gz` (or similar)
   - Example: `Sample1_R1.fastq.gz` and `Sample1_R2.fastq.gz`

2. **In Design Matrices and Contrast Tables**:
   - Include only the sample name (without the `_R1`/`_R2` suffix)
   - The pipeline will automatically find and process both read files together
   - Multiple replicates of the same condition should be listed separately in the design matrix

### Example Design Matrix for Paired-End Data

```
Sample,Treatment,Control
Control_Rep1,0,1
Control_Rep2,0,1
Treatment_Rep1,1,0
Treatment_Rep2,1,0
```

### Example Contrast Table for Paired-End Data

```
contrast,treatment,control
experiment1,Treatment_Rep1,Treatment_Rep2,Control_Rep1,Control_Rep2
```

Note that in the contrast table, multiple replicates of the same condition are comma-separated in the treatment or control column.