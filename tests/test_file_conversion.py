#!/usr/bin/env python3
"""
Unit tests for the file conversion utility.
"""

import os
import sys
import tempfile
import unittest
import pandas as pd
from pathlib import Path

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

try:
    from core.file_handling import convert_file_to_tab_delimited
    from convert_input_files import find_input_files, convert_files
except ImportError:
    # Try with analysis_pipeline prefix
    from analysis_pipeline.core.file_handling import convert_file_to_tab_delimited
    from analysis_pipeline.convert_input_files import find_input_files, convert_files


class TestFileConversion(unittest.TestCase):
    """Test cases for file conversion utility."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory
        self.test_dir = tempfile.TemporaryDirectory()
        self.temp_dir_path = Path(self.test_dir.name)
        
        # Create test design matrix CSV
        self.design_matrix_path = self.temp_dir_path / "design_matrix.csv"
        design_matrix_data = pd.DataFrame({
            'Sample': ['sample1', 'sample2', 'sample3', 'sample4'],
            'Control': [1, 1, 0, 0],
            'Treatment': [0, 0, 1, 1]
        })
        design_matrix_data.to_csv(self.design_matrix_path, index=False)
        
        # Create test contrast table CSV (standard format)
        self.contrast_table_path = self.temp_dir_path / "contrast_table.csv"
        contrast_table_data = pd.DataFrame({
            'contrast': ['test_contrast', 'test_contrast2'],
            'control': ['sample1,sample2', 'sample1'],
            'treatment': ['sample3,sample4', 'sample3']
        })
        contrast_table_data.to_csv(self.contrast_table_path, index=False)
        
        # Create test contrast table CSV with duplicate columns
        self.contrast_duplicate_cols_path = self.temp_dir_path / "contrast_duplicate_cols.csv"
        with open(self.contrast_duplicate_cols_path, 'w') as f:
            f.write("contrast,control,control,treatment,treatment\n")
            f.write("test_contrast,control1,control2,treatment1,treatment2\n")
            f.write("test_contrast2,control3,control4,treatment3,treatment4\n")
        
        # Create test contrast table with whitespace issues
        self.contrast_whitespace_path = self.temp_dir_path / "contrast_whitespace.csv"
        with open(self.contrast_whitespace_path, 'w') as f:
            f.write("contrast,control,treatment\n")
            f.write("test_contrast, sample1 , sample3 \n")
            f.write("test_contrast2,  sample1,sample2  ,  sample3,sample4  \n")
        
        # Create test contrast table with whitespace issues and duplicate columns
        self.contrast_whitespace_duplicate_path = self.temp_dir_path / "contrast_whitespace_duplicate.csv"
        with open(self.contrast_whitespace_duplicate_path, 'w') as f:
            f.write("contrast,control,control,treatment,treatment\n")
            f.write("test_contrast, control1 , control2 , treatment1 , treatment2 \n")
            f.write("test_contrast2,  control3,  control4,  treatment3,  treatment4\n")
        
        # Create another test file that doesn't match patterns
        self.other_file_path = self.temp_dir_path / "other_file.csv"
        other_data = pd.DataFrame({'A': [1, 2], 'B': [3, 4]})
        other_data.to_csv(self.other_file_path, index=False)

        # Set up experiment-specific directory structure
        self.exp1_dir = self.temp_dir_path / "experiment1"
        self.exp1_dir.mkdir(exist_ok=True)
        
        self.exp2_dir = self.temp_dir_path / "experiment2"
        self.exp2_dir.mkdir(exist_ok=True)
        
        # Create files in experiment1 directory
        self.exp1_design_matrix_path = self.exp1_dir / "design_matrix.csv"
        design_matrix_data.to_csv(self.exp1_design_matrix_path, index=False)
        
        self.exp1_contrast_table_path = self.exp1_dir / "contrasts.csv"
        contrast_table_data.to_csv(self.exp1_contrast_table_path, index=False)
        
        # Create files in experiment2 directory
        self.exp2_design_matrix_path = self.exp2_dir / "design_matrix.csv"
        design_matrix_data.to_csv(self.exp2_design_matrix_path, index=False)
        
        self.exp2_contrast_table_path = self.exp2_dir / "contrast_table.csv"
        contrast_table_data.to_csv(self.exp2_contrast_table_path, index=False)

    def tearDown(self):
        """Tear down test fixtures."""
        self.test_dir.cleanup()

    def test_find_input_files(self):
        """Test finding input files in a directory."""
        files = find_input_files(self.temp_dir_path)
        
        # Should find design matrices and contrast tables
        self.assertGreaterEqual(len(files['design_matrices']), 1)
        self.assertGreaterEqual(len(files['contrast_tables']), 2)
        
        # Check file paths
        self.assertIn(str(self.design_matrix_path), files['design_matrices'])
        self.assertIn(str(self.contrast_table_path), files['contrast_tables'])
        self.assertIn(str(self.contrast_duplicate_cols_path), files['contrast_tables'])

    def test_find_input_files_with_experiment(self):
        """Test finding input files in a specific experiment directory."""
        # Test with experiment1
        files = find_input_files(self.temp_dir_path, "experiment1")
        
        # Should find only experiment1's files
        self.assertEqual(len(files['design_matrices']), 1)
        self.assertEqual(len(files['contrast_tables']), 1)
        
        # Check file paths
        self.assertIn(str(self.exp1_design_matrix_path), files['design_matrices'])
        self.assertIn(str(self.exp1_contrast_table_path), files['contrast_tables'])
        
        # Test with experiment2
        files = find_input_files(self.temp_dir_path, "experiment2")
        
        # Should find only experiment2's files
        self.assertEqual(len(files['design_matrices']), 1)
        self.assertEqual(len(files['contrast_tables']), 1)
        
        # Check file paths
        self.assertIn(str(self.exp2_design_matrix_path), files['design_matrices'])
        self.assertIn(str(self.exp2_contrast_table_path), files['contrast_tables'])

    def test_find_input_files_with_nonexistent_experiment(self):
        """Test finding input files with a nonexistent experiment name."""
        files = find_input_files(self.temp_dir_path, "nonexistent_experiment")
        
        # Should fall back to base directory
        self.assertGreaterEqual(len(files['design_matrices']), 1)
        self.assertGreaterEqual(len(files['contrast_tables']), 2)
        
        # Check that it includes the base directory files
        self.assertIn(str(self.design_matrix_path), files['design_matrices'])

    def test_convert_design_matrix(self):
        """Test converting a design matrix file."""
        # Convert the design matrix
        converted_file = convert_file_to_tab_delimited(str(self.design_matrix_path))
        
        # Converted file should exist
        self.assertTrue(os.path.exists(converted_file))
        
        # Check that the file has .txt extension
        self.assertTrue(converted_file.endswith('.txt'))
        
        # Read the converted file and check the format
        df = pd.read_csv(converted_file, sep='\t')
        
        # Check that all expected columns are present
        self.assertIn('Sample', df.columns)
        self.assertIn('Control', df.columns)
        self.assertIn('Treatment', df.columns)
        
        # Check number of rows
        self.assertEqual(len(df), 4)

    def test_convert_contrast_table(self):
        """Test converting a contrast table file."""
        # Convert the contrast table
        converted_file = convert_file_to_tab_delimited(str(self.contrast_table_path))
        
        # Converted file should exist
        self.assertTrue(os.path.exists(converted_file))
        
        # Check that the file has .txt extension
        self.assertTrue(converted_file.endswith('.txt'))
        
        # Read the converted file and check the format
        df = pd.read_csv(converted_file, sep='\t')
        
        # Check that all required columns are present
        self.assertIn('contrast', df.columns)
        self.assertIn('control', df.columns)
        self.assertIn('treatment', df.columns)
        
        # Check that the data is correct
        self.assertEqual(df.iloc[0]['contrast'], 'test_contrast')
        self.assertEqual(df.iloc[0]['control'], 'sample1,sample2')
        self.assertEqual(df.iloc[0]['treatment'], 'sample3,sample4')

    def test_convert_contrast_table_with_duplicate_columns(self):
        """Test converting a contrast table file with duplicate column names."""
        # Convert the contrast table with duplicate columns
        converted_file = convert_file_to_tab_delimited(str(self.contrast_duplicate_cols_path))
        
        # Converted file should exist
        self.assertTrue(os.path.exists(converted_file))
        
        # Check that the file has .txt extension
        self.assertTrue(converted_file.endswith('.txt'))
        
        # Read the converted file and check the format
        df = pd.read_csv(converted_file, sep='\t')
        
        # Check that all required columns are present and only appear once
        self.assertIn('contrast', df.columns)
        self.assertIn('control', df.columns)
        self.assertIn('treatment', df.columns)
        self.assertEqual(len(df.columns), 3)  # Only 3 unique columns
        
        # Check that the duplicate columns were properly joined
        self.assertEqual(df.iloc[0]['contrast'], 'test_contrast')
        self.assertEqual(df.iloc[0]['control'], 'control1,control2')
        self.assertEqual(df.iloc[0]['treatment'], 'treatment1,treatment2')
        
        self.assertEqual(df.iloc[1]['contrast'], 'test_contrast2')
        self.assertEqual(df.iloc[1]['control'], 'control3,control4')
        self.assertEqual(df.iloc[1]['treatment'], 'treatment3,treatment4')

    def test_convert_files_function(self):
        """Test the convert_files function that handles multiple files."""
        # Convert both files
        files_to_convert = [str(self.design_matrix_path), str(self.contrast_table_path)]
        converted_files = convert_files(files_to_convert)
        
        # Should have converted both files
        self.assertEqual(len(converted_files), 2)
        
        # Check that all converted files exist
        for file_path in converted_files:
            self.assertTrue(os.path.exists(file_path))
            self.assertTrue(file_path.endswith('.txt'))

    def test_convert_to_custom_output_dir(self):
        """Test converting files to a custom output directory."""
        # Create a custom output directory
        output_dir = self.temp_dir_path / "output"
        output_dir.mkdir(exist_ok=True)
        
        # Convert a file to the custom output directory
        converted_file = convert_file_to_tab_delimited(
            str(self.design_matrix_path), 
            output_dir=str(output_dir)
        )
        
        # Check that the file was saved in the custom output directory
        self.assertTrue(os.path.exists(converted_file))
        self.assertTrue(str(output_dir) in converted_file)

    def test_convert_contrast_table_with_whitespace(self):
        """Test converting a contrast table file with whitespace issues."""
        # Convert the contrast table with whitespace issues
        converted_file = convert_file_to_tab_delimited(str(self.contrast_whitespace_path))
        
        # Converted file should exist
        self.assertTrue(os.path.exists(converted_file))
        
        # Read the converted file and check the format
        df = pd.read_csv(converted_file, sep='\t')
        
        # Check that whitespace was cleaned from all values
        self.assertEqual(df.iloc[0]['control'], 'sample1')
        self.assertEqual(df.iloc[0]['treatment'], 'sample3')
        
        # Check that whitespace was cleaned from comma-separated lists
        self.assertEqual(df.iloc[1]['control'], 'sample1,sample2')
        self.assertEqual(df.iloc[1]['treatment'], 'sample3,sample4')

    def test_convert_contrast_table_with_whitespace_and_duplicate_columns(self):
        """Test converting a contrast table with both whitespace and duplicate columns."""
        # Convert the contrast table with whitespace and duplicate columns
        converted_file = convert_file_to_tab_delimited(str(self.contrast_whitespace_duplicate_path))
        
        # Converted file should exist
        self.assertTrue(os.path.exists(converted_file))
        
        # Read the converted file and check the format
        df = pd.read_csv(converted_file, sep='\t')
        
        # Check that whitespace was cleaned and duplicate columns were properly joined
        self.assertEqual(df.iloc[0]['control'], 'control1,control2')
        self.assertEqual(df.iloc[0]['treatment'], 'treatment1,treatment2')
        self.assertEqual(df.iloc[1]['control'], 'control3,control4')
        self.assertEqual(df.iloc[1]['treatment'], 'treatment3,treatment4')

    def test_convert_files_in_experiment_dir(self):
        """Test converting files in an experiment-specific directory."""
        # Find files in experiment1
        files = find_input_files(self.temp_dir_path, "experiment1")
        
        # Convert all found files
        converted_files = convert_files(
            files['design_matrices'] + files['contrast_tables']
        )
        
        # Should have converted both files from experiment1
        self.assertEqual(len(converted_files), 2)
        
        # Ensure the converted files exist and have .txt extension
        for file_path in converted_files:
            self.assertTrue(os.path.exists(file_path))
            self.assertTrue(file_path.endswith('.txt'))
        
        # Verify that one of the files is the design matrix
        design_matrix_converted = any(
            os.path.basename(f).startswith("design_matrix") for f in converted_files
        )
        self.assertTrue(design_matrix_converted)
        
        # Verify that one of the files is the contrast table
        contrast_table_converted = any(
            os.path.basename(f).startswith("contrasts") for f in converted_files
        )
        self.assertTrue(contrast_table_converted)


if __name__ == '__main__':
    unittest.main() 