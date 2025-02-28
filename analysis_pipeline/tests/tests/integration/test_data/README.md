# Integration Test Data

This directory is used for storing test data files for integration testing. 

The test data files are created dynamically by the test scripts in the `test_runners` directory. There is no need to manually add files here, as the tests will create the necessary files when they are run.

## Test Data Types

The integration tests will create several types of test data:

1. **FASTQ Files**: Small, synthetic FASTQ files containing minimal read data for testing.
2. **Library Files**: sgRNA library files with gene mappings.
3. **Contrasts Files**: Files defining experimental contrasts.
4. **Empty Files**: For testing error handling of empty input.
5. **Malformed Files**: For testing error handling of malformed input.

## Directory Structure

The directory structure created during testing will typically look like:

```
test_data/
  ├── sample1.fastq
  ├── sample2.fastq 
  ├── empty_sample.fastq
  ├── malformed_sample.fastq
  ├── test_library.txt
  └── test_contrasts.txt
```

## Note

This directory and its contents are meant for testing only. The files are simplified versions of real data files and are not suitable for actual analysis. 