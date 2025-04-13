# Input File Conversion Utility

This utility helps convert CSV-formatted design matrices and contrast tables to the tab-delimited format required by MAGeCK and DrugZ analysis tools.

## Features

- Converts CSV files to tab-delimited text files
- Automatically detects design matrices and contrast tables
- Handles different contrast table formats (standard and multiple column)
- Cleans unwanted whitespace from values to prevent validation errors
- Preserves data structure while ensuring compatibility with analysis tools
- **NEW: Supports experiment-specific directory structures**
- **NEW: Integrated with Snakemake workflow**

## Requirements

- Python 3.6 or higher
- pandas

## Usage

### Converting a Single File

```bash
python convert_input_files.py --file /path/to/your/file.csv
```

### Converting All Files in a Directory

```bash
python convert_input_files.py --dir /path/to/your/directory
```

Or simply:

```bash
python convert_input_files.py /path/to/your/directory
```

### Specifying an Experiment

For multi-experiment directory structures, you can specify which experiment's files to convert:

```bash
python convert_input_files.py --dir /path/to/screens_directory --experiment experiment1
```

This will search for files specifically in `/path/to/screens_directory/experiment1/`.

### Specifying an Output Directory

By default, converted files are saved in the same directory as the input files. You can specify a different output directory:

```bash
python convert_input_files.py --file /path/to/your/file.csv --output /path/to/output/directory
```

### Additional Options

- `--verbose` or `-v`: Enable verbose logging
- `--help` or `-h`: Display help information

## Directory Structure Support

The utility is designed to work with both the new experimental directory structure and legacy formats:

### New Experiment-Specific Structure

```
<screens_directory>/
├── <experiment1>/           # First experiment directory
│   ├── library.csv          # CRISPR library file
│   ├── contrasts.csv        # Contrast definitions file
│   └── design_matrix.csv    # Design matrix for MLE analysis
│
├── <experiment2>/           # Second experiment directory
│   ├── library.csv
│   ├── contrasts.csv
│   └── design_matrix.csv
```

When run with the `--experiment` flag, the utility will search specifically in that experiment's subdirectory. When run without this flag, it will search all experiment directories.

### Legacy Structure

```
<input_dir>/
├── library.csv              # CRISPR library file
├── contrasts.csv            # Contrast definitions file 
├── design_matrix.csv        # Design matrix for MLE analysis
```

The utility will automatically detect and process files in legacy directory structures.

## Integration with Snakemake

This utility is integrated with the Snakemake workflow in the CRISPR Analysis Pipeline. The conversion happens automatically as part of the workflow:

1. Snakemake detects CSV files that need to be converted
2. The conversion utility transforms these files to tab-delimited format
3. Analysis proceeds with the converted files

You don't need to run the conversion utility manually when using the Snakemake workflow - it happens automatically. However, you may want to run it separately for pre-processing your input files before analysis.

## File Format Requirements

### Design Matrix Files

The script looks for CSV files with names containing "design" or "design_matrix". These files should have:

- A 'Sample' column with sample names
- One or more condition columns
- CSV format (comma-separated)

Example design matrix CSV format:
```
Sample,Condition1,Condition2
sample1,1,0
sample2,0,1
```

The script will convert this to a tab-separated text file:
```
Sample	Condition1	Condition2
sample1	1	0
sample2	0	1
```

### Contrast Table Files

The script supports two formats for contrast table CSV files:

#### Standard Format

The standard format has one column each for 'contrast', 'control', and 'treatment':

```
contrast,control,treatment
test_contrast,control1,treatment1
test_contrast2,control1,treatment1,treatment2
test_contrast3,control1,control2,treatment1
```

Where multiple samples in a single cell are comma-separated.

#### Multiple Column Format

For easier human editing, you can also use a format with multiple 'control' or 'treatment' columns:

```
contrast,control,control,treatment,treatment
test_contrast,control1,control2,treatment1,treatment2
```

The utility will automatically detect duplicate columns and merge them, creating a properly formatted tab-delimited file with comma-separated sample lists:

```
contrast	control	treatment
test_contrast	control1,control2	treatment1,treatment2
```

This is particularly useful when:
- You have many samples and want to keep them visually organized
- You're creating files in a spreadsheet application where comma-separated values in a single cell might be confusing

The script will convert both formats to the tab-separated text format required by MAGeCK and DrugZ:

```
contrast	control	treatment
test_contrast	control1,control2	treatment1,treatment2
```

## Additional Processing

The utility performs several helpful transformations on your input files:

### Whitespace Cleaning

The tool automatically cleans unwanted whitespace from all values. This prevents validation errors that can occur when sample names have extra spaces. For example:

Input (with unwanted spaces):
```
contrast,control,treatment
test_contrast, sample1 , sample3 
test_contrast2,  sample1,sample2  ,  sample3,sample4  
```

Output (cleaned):
```
contrast	control	treatment
test_contrast	sample1	sample3
test_contrast2	sample1,sample2	sample3,sample4
```

This is especially important for sample names, as the analysis pipeline needs to match exact sample names when validating files.

## Output

All converted files will be saved with a `.txt` extension. The script will create this file in:
- The same directory as the input file (default)
- A specified output directory (if provided)

The script logs all actions to the console, including:
- Number of files found
- Conversion successes and failures
- Final conversion summary

## Pipeline Input Requirements

The CRISPR Analysis Pipeline handles different types of input data, which affects the required files:

### When Starting with FASTQ Files

If you have a `fastq` subdirectory in your experiment directory:

- **Library file is REQUIRED** in the experiment directory
- The pipeline needs the library file to map sgRNA sequences to genes
- The utility will convert your design matrix and contrast table files

Example directory structure:
```
experiment_dir/
  ├── fastq/
  │   ├── Sample1_R1.fastq.gz
  │   ├── Sample1_R2.fastq.gz  # For paired-end data
  │   └── ...
  ├── library.csv          # REQUIRED for FASTQ processing
  ├── design_matrix.csv    # Will be converted to .txt
  └── contrast_table.csv   # Will be converted to .txt
```

### When Starting with Count Files

If you have a `counts` subdirectory in your experiment directory:

- **Library file is NOT required** as the sgRNA-to-gene mapping has already been done
- The utility will still convert your design matrix and contrast table files

Example directory structure:
```
experiment_dir/
  ├── counts/
  │   ├── count_table.csv
  │   └── ...
  ├── design_matrix.csv    # Will be converted to .txt
  └── contrast_table.csv   # Will be converted to .txt
```

### Mixed Input Directories

For directories containing both FASTQ and count files:

- **Library file is required** as the FASTQ files need to be processed
- The pipeline will use the appropriate input type for each sample

Example directory structure:
```
experiment_dir/
  ├── fastq/
  │   ├── Sample1_R1.fastq.gz
  │   └── ...
  ├── counts/
  │   ├── Sample2.count
  │   └── ...
  ├── library.csv        # REQUIRED for processing FASTQ files
  ├── design_matrix.csv  # Will be converted to .txt
  └── contrast_table.csv # Will be converted to .txt
``` 