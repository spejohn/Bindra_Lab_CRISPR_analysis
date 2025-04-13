# CRISPR Analysis Pipeline

A modular pipeline for CRISPR screening data analysis using Snakemake, supporting MAGeCK (RRA and MLE methods) and DrugZ analysis.

**Current Version: 0.1.0**

## Overview

This pipeline provides a robust and reproducible workflow for analyzing CRISPR screen data, starting from either raw FASTQ files or pre-computed count tables. It automates quality control, read counting (if necessary), normalization, and statistical analysis using established tools like MAGeCK and DrugZ within Docker containers.

## Features

-   **Snakemake Workflow:** Ensures reproducibility and efficient execution.
-   **MAGeCK & DrugZ Integration:** Supports standard analysis methods.
-   **FASTQ or Count Input:** Handles both raw sequencing data and count tables.
-   **Dockerized Tools:** Encapsulates dependencies for consistent execution.
-   **Parallel Processing:** Leverages multiple cores for faster analysis.
-   **Quality Control:** Includes steps for assessing data quality.
-   **Standardized Outputs:** Generates organized results and reports.
-   **Cross-Platform:** Runs on Windows, macOS, and Linux.

## Quick Start

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/your-username/crispr_analysis_pipeline.git
    cd crispr_analysis_pipeline
    ```

2.  **Set Up Environment (Recommended):**
    ```bash
    python -m venv venv
    # Activate:
    # Windows: .\venv\Scripts\activate
    # macOS/Linux: source venv/bin/activate
    ```

3.  **Install Dependencies:**
    ```bash
    pip install . # Or pip install -e . 
    ```
    *Ensure Docker is installed and running.*

4.  **Prepare Input Data:**
    Organize your input files (FASTQ/counts, library, contrasts) according to the structure described in [Input and Output Structure](docs/INPUT_OUTPUT.md).
    *Use the [File Conversion Utility](docs/FILE_CONVERSION.md) if needed to convert CSV contrast/design files to TXT.*

5.  **Run the Pipeline (Snakemake Runner):**
    ```bash
    python run_snakemake.py /path/to/your/input_directory -o /path/to/output -j <num_cores>
    ```
    The runner automatically detects experiments and input types.

## Documentation

For more detailed information, please refer to the documentation in the `docs/` directory:

-   **[Installation Guide](docs/INSTALL.md):** Detailed setup instructions, including development mode and Docker configuration.
-   **[Input and Output Structure](docs/INPUT_OUTPUT.md):** Comprehensive guide on required input files, formats, directory layouts, and output organization.
-   **[Testing Guide](docs/TESTING.md):** Instructions for running various tests to verify the pipeline.
-   **[File Conversion Utility](docs/FILE_CONVERSION.md):** How to use the utility for converting input CSV files.
-   **[Downloading Data Guide](docs/DOWNLOADING_DATA.md):** Information related to downloading data from public repositories (e.g., SRA).

## Main Scripts

-   `run_snakemake.py`: Recommended entry point for running the Snakemake workflow.
-   `run_multiple_screens.py`: (Potentially for running multiple experiments outside Snakemake - verify usage)
-   `run_qc_analysis.py`: (Potentially for running QC standalone - verify usage)
-   `run_integration_tests.py`: Executes integration tests.

Utility scripts can be found in the `scripts/` directory.

## Contributing

Contributions are welcome! Please refer to the [Installation Guide](docs/INSTALL.md) for development setup and the [Testing Guide](docs/TESTING.md) for verifying changes.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.

## License

(Add License Information Here)
