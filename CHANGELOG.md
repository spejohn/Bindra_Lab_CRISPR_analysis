# Changelog

All notable changes to the CRISPR Analysis Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-04-13

### Changed
- Refactored project structure significantly for clarity and maintainability.
- Consolidated all documentation (.md, .txt files) into a root `README.md` and a dedicated `docs/` directory.
- Moved utility scripts (`install_dev.py`, `download_fastq_files.py`, `profile_workflow.py`) into a dedicated `scripts/` directory.
- Moved test scripts (`test_*.py`) into the `tests/` directory.
- Removed redundant `analysis_pipeline/` directory containing duplicate code, configuration, and documentation.
- Corrected import paths in `analysis/mageck_analysis.py` to reflect the new package structure.
- Updated `setup.py` to correctly find package modules (`core`, `analysis`, `docker`, etc.) in the root and removed incorrect/obsolete entry points.
- Updated `.gitignore` to explicitly ignore `cleanup_backups/` directory.

### Removed
- Deleted redundant `cleanup_backups/` directory.
- Deleted redundant `analysis_pipeline/Snakefile` and `analysis_pipeline/config.yaml`.
- Removed `convert_input_files` script entry point from `setup.py` (functionality remains as internal utility).

## [0.1.0] - 2023-11-15

### Added
- Initial stable version with standardized function names
- Version management via `__version__` constant
- Support for parallel processing in the counting step
- Automatic detection of Snakemake environment to disable internal parallelism
- Progress reporting throughout the pipeline execution

### Changed
- Removed backward compatibility wrapper functions
- Standardized function interfaces and parameter names
- Improved consistency across the codebase
- Enhanced error handling with retry mechanisms
- Updated documentation and README

### Fixed
- Fixed failing benchmark tests
- Addressed issues with path handling in different operating systems
- Resolved race conditions in parallel code execution 