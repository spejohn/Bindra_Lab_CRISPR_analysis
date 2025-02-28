# Changelog

All notable changes to the CRISPR Analysis Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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