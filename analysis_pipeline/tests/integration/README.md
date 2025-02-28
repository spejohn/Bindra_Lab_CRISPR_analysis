# CRISPR Analysis Pipeline - Integration Tests

This directory contains integration tests for the CRISPR analysis pipeline using Snakemake to orchestrate the tests.

## Overview

The integration tests use Snakemake to run complete workflows with test data, verifying that:

1. All steps of the pipeline run successfully with different configuration options
2. Output files are generated correctly and contain expected data
3. The pipeline correctly detects new input files and processes them
4. Error handling works as expected with problematic input files

## Test Structure

- `test_workflow/` - Contains test-specific Snakefile and configuration
- `test_data/` - Minimal test datasets for running integration tests
- `test_runners/` - Python scripts for launching and verifying the tests
- `conftest.py` - Pytest fixtures for integration testing

## Running Tests

To run all integration tests:

```bash
pytest tests/integration
```

To run a specific test workflow:

```bash
python tests/integration/test_runners/test_full_workflow.py
```

## Test Data

The test data is minimal to allow fast test execution but covers different scenarios:
- Normal FASTQ files with expected formats
- Empty or malformed files to test error handling
- Multiple screens with different layouts to test organization
- Various library files to test compatibility

## Test Workflows

The tests cover several workflows:
1. Complete pipeline execution (end-to-end test)
2. Incremental workflow testing (adding new samples to existing results)
3. Error condition testing (malformed inputs, missing files)
4. Configuration variations (different normalization methods, Docker usage) 