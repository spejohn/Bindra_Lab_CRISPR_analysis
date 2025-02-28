#!/usr/bin/env python3
"""
Unit tests for the benchmarking pipeline.
"""

import os
import pandas as pd
import pytest
from unittest.mock import patch, mock_open, MagicMock
from pathlib import Path

from analysis_pipeline.benchmarking.benchmark_pipeline import (
    find_file_pairs,
    copy_files_for_qc,
    run_direct_comparison
)


@pytest.fixture
def mock_directory_structure():
    """Create a mock directory structure for testing."""
    return {
        "generated_dir": {
            "experiment1": {
                "contrast1": {
                    "screen1_gDZ.csv": "GENE,normZ,FDR\ngene1,2.5,0.01\ngene2,-1.8,0.05"
                },
                "contrast2": {
                    "screen2_gDZ.csv": "GENE,normZ,FDR\ngene1,1.9,0.02\ngene2,-2.1,0.03"
                }
            },
            "experiment2": {
                "contrast1": {
                    "screen1_gMGK.csv": "id,neg|lfc,neg|p-value,neg|fdr\ngene1,-0.8,0.01,0.03\ngene2,1.2,0.02,0.04"
                }
            }
        },
        "published_dir": {
            "experiment1": {
                "screen1_pDZ.csv": "GENE,normZ,FDR\ngene1,2.4,0.012\ngene2,-1.75,0.055",
                "screen2_pDZ.csv": "GENE,normZ,FDR\ngene1,1.85,0.022\ngene2,-2.05,0.033",
            },
            "experiment2": {
                "screen1_pMGK.csv": "id,neg|lfc,neg|p-value,neg|fdr\ngene1,-0.78,0.012,0.032\ngene2,1.18,0.022,0.042"
            }
        }
    }


@patch('analysis_pipeline.benchmarking.benchmark_pipeline.Path')
def test_find_file_pairs(mock_path, mock_directory_structure):
    """Test finding matching file pairs between directories."""
    # Setup mock paths and directory structure
    mock_path_instance = MagicMock()
    mock_gen_dir = MagicMock()
    mock_pub_dir = MagicMock()
    
    # Setup mock iterdir for experiments
    mock_gen_exp1 = MagicMock()
    mock_gen_exp1.name = "experiment1"
    mock_gen_exp1.is_dir.return_value = True
    mock_gen_exp1.exists.return_value = True
    
    mock_gen_exp2 = MagicMock()
    mock_gen_exp2.name = "experiment2"
    mock_gen_exp2.is_dir.return_value = True
    mock_gen_exp2.exists.return_value = True
    
    mock_gen_dir.iterdir.return_value = [mock_gen_exp1, mock_gen_exp2]
    
    # Setup mock experiments in published dir
    mock_pub_exp1 = MagicMock()
    mock_pub_exp1.name = "experiment1"
    mock_pub_exp1.exists.return_value = True
    mock_pub_exp1.is_dir.return_value = True
    
    mock_pub_exp2 = MagicMock()
    mock_pub_exp2.name = "experiment2"
    mock_pub_exp2.exists.return_value = True
    mock_pub_exp2.is_dir.return_value = True
    
    # Setup paths for contrasts in generated dir
    mock_gen_exp1_cont1 = MagicMock()
    mock_gen_exp1_cont1.name = "contrast1"
    mock_gen_exp1_cont1.is_dir.return_value = True
    mock_gen_exp1_cont1.exists.return_value = True
    
    mock_gen_exp1_cont2 = MagicMock()
    mock_gen_exp1_cont2.name = "contrast2"
    mock_gen_exp1_cont2.is_dir.return_value = True
    mock_gen_exp1_cont2.exists.return_value = True
    
    mock_gen_exp1.iterdir.return_value = [mock_gen_exp1_cont1, mock_gen_exp1_cont2]
    
    mock_gen_exp2_cont1 = MagicMock()
    mock_gen_exp2_cont1.name = "contrast1"
    mock_gen_exp2_cont1.is_dir.return_value = True
    mock_gen_exp2_cont1.exists.return_value = True
    
    mock_gen_exp2.iterdir.return_value = [mock_gen_exp2_cont1]
    
    # Setup files in each contrast directory
    mock_screen1_gdz = MagicMock()
    mock_screen1_gdz.name = "screen1_gDZ.csv"
    
    mock_screen2_gdz = MagicMock()
    mock_screen2_gdz.name = "screen2_gDZ.csv"
    
    mock_screen1_gmgk = MagicMock()
    mock_screen1_gmgk.name = "screen1_gMGK.csv"
    
    mock_gen_exp1_cont1.glob.side_effect = lambda pattern: [mock_screen1_gdz] if "_gDZ.csv" in pattern else []
    mock_gen_exp1_cont2.glob.side_effect = lambda pattern: [mock_screen2_gdz] if "_gDZ.csv" in pattern else []
    mock_gen_exp2_cont1.glob.side_effect = lambda pattern: [mock_screen1_gmgk] if "_gMGK.csv" in pattern else []
    
    # Setup published files and their paths
    mock_pub_exp1_cont1 = MagicMock()
    mock_pub_exp1_cont1.name = "contrast1"
    mock_pub_exp1_cont1.exists.return_value = True
    mock_pub_exp1_cont1.is_dir.return_value = True
    
    mock_pub_exp1_cont2 = MagicMock()
    mock_pub_exp1_cont2.name = "contrast2"
    mock_pub_exp1_cont2.exists.return_value = True
    mock_pub_exp1_cont2.is_dir.return_value = True
    
    mock_pub_exp1.iterdir.return_value = [mock_pub_exp1_cont1, mock_pub_exp1_cont2]
    
    mock_pub_exp2_cont1 = MagicMock()
    mock_pub_exp2_cont1.name = "contrast1"
    mock_pub_exp2_cont1.exists.return_value = True
    mock_pub_exp2_cont1.is_dir.return_value = True
    
    mock_pub_exp2.iterdir.return_value = [mock_pub_exp2_cont1]
    
    # Setup published files
    mock_screen1_pdz = MagicMock()
    mock_screen1_pdz.name = "screen1_pDZ.csv"
    
    mock_screen2_pdz = MagicMock()
    mock_screen2_pdz.name = "screen2_pDZ.csv"
    
    mock_screen1_pmgk = MagicMock()
    mock_screen1_pmgk.name = "screen1_pMGK.csv"
    
    # Mock the exists check for the published files
    mock_pub_exist_true = MagicMock()
    mock_pub_exist_true.exists.return_value = True
    
    # Setup Path return values
    def mock_path_call(arg):
        if arg == "generated_dir":
            return mock_gen_dir
        elif arg == "published_dir":
            return mock_pub_dir
        elif arg == mock_gen_dir / "experiment1":
            return mock_gen_exp1
        elif arg == mock_gen_dir / "experiment2":
            return mock_gen_exp2
        elif arg == mock_pub_dir / "experiment1":
            return mock_pub_exp1
        elif arg == mock_pub_dir / "experiment2":
            return mock_pub_exp2
        elif arg == mock_gen_exp1 / "contrast1":
            return mock_gen_exp1_cont1
        elif arg == mock_gen_exp1 / "contrast2":
            return mock_gen_exp1_cont2
        elif arg == mock_gen_exp2 / "contrast1":
            return mock_gen_exp2_cont1
        elif arg == mock_pub_exp1 / "contrast1":
            return mock_pub_exp1_cont1
        elif arg == mock_pub_exp1 / "contrast2":
            return mock_pub_exp1_cont2
        elif arg == mock_pub_exp2 / "contrast1":
            return mock_pub_exp2_cont1
        # Handle published file paths in both experiment and contrast directories
        elif arg == mock_pub_exp1 / "screen1_pDZ.csv":
            result = MagicMock()
            result.name = "screen1_pDZ.csv"
            result.exists.return_value = True
            return result
        elif arg == mock_pub_exp1 / "screen2_pDZ.csv":
            result = MagicMock()
            result.name = "screen2_pDZ.csv"
            result.exists.return_value = True
            return result
        elif arg == mock_pub_exp2 / "screen1_pMGK.csv":
            result = MagicMock()
            result.name = "screen1_pMGK.csv"
            result.exists.return_value = True
            return result
        elif arg == mock_pub_exp1 / "contrast1" / "screen1_pDZ.csv":
            result = MagicMock()
            result.name = "screen1_pDZ.csv"
            result.exists.return_value = False  # Not found in contrast dir, will fall back to exp dir
            return result
        elif arg == mock_pub_exp1 / "contrast2" / "screen2_pDZ.csv":
            result = MagicMock()
            result.name = "screen2_pDZ.csv" 
            result.exists.return_value = False  # Not found in contrast dir, will fall back to exp dir
            return result
        elif arg == mock_pub_exp2 / "contrast1" / "screen1_pMGK.csv":
            result = MagicMock()
            result.name = "screen1_pMGK.csv"
            result.exists.return_value = False  # Not found in contrast dir, will fall back to exp dir
            return result
        return MagicMock()
    
    mock_path.side_effect = mock_path_call
    
    # Call the function
    result = find_file_pairs(
        generated_dir="generated_dir",
        published_dir="published_dir"
    )
    
    # Check results
    assert len(result) == 3
    
    # Check DrugZ file pairs
    drugz_pairs = [pair for pair in result if pair["analysis_type"] == "DrugZ"]
    assert len(drugz_pairs) == 2
    
    # Check MAGeCK file pairs
    mageck_pairs = [pair for pair in result if pair["analysis_type"] == "MAGeCK"]
    assert len(mageck_pairs) == 1


@patch('os.makedirs')
@patch('shutil.copy2')
@patch('analysis_pipeline.benchmarking.benchmark_pipeline.ensure_output_dir')
@patch('analysis_pipeline.benchmarking.benchmark_pipeline.Path')
def test_copy_files_for_qc(mock_path, mock_ensure_dir, mock_copy, mock_makedirs):
    """Test copying files to the QC directory structure."""
    # Setup file pairs
    file_pairs = [
        {
            "experiment": "experiment1",
            "contrast": "contrast1",
            "analysis_type": "DrugZ",
            "files": {
                "generated": "/path/to/generated/experiment1/contrast1/screen1_gDZ.csv",
                "published": "/path/to/published/experiment1/screen1_pDZ.csv"
            }
        },
        {
            "experiment": "experiment2",
            "contrast": "contrast1",
            "analysis_type": "MAGeCK",
            "files": {
                "generated": "/path/to/generated/experiment2/contrast1/screen1_gMGK.csv",
                "published": "/path/to/published/experiment2/screen1_pMGK.csv"
            }
        }
    ]
    
    # Mock Path objects
    mock_path.side_effect = lambda arg: MagicMock(name=str(arg))
    
    # Call function
    result = copy_files_for_qc(file_pairs, "/output/dir")
    
    # Check results
    assert result["DrugZ"]["generated"] == 1
    assert result["DrugZ"]["published"] == 1
    assert result["MAGeCK"]["generated"] == 1
    assert result["MAGeCK"]["published"] == 1
    
    # Verify copies were made
    assert mock_copy.call_count == 4


@patch('analysis_pipeline.benchmarking.benchmark_pipeline.ensure_output_dir')
@patch('analysis_pipeline.benchmarking.benchmark_pipeline.calc_stats')
@patch('pandas.read_csv')
@patch('pandas.DataFrame.to_csv')
def test_run_direct_comparison(mock_to_csv, mock_read_csv, mock_calc_stats, mock_ensure_dir):
    """Test running a direct comparison between file pairs."""
    # Setup file pairs
    file_pairs = [
        {
            "experiment": "experiment1",
            "contrast": "contrast1",
            "analysis_type": "DrugZ",
            "files": {
                "generated": "/path/to/generated/experiment1/contrast1/screen1_gDZ.csv",
                "published": "/path/to/published/experiment1/screen1_pDZ.csv"
            }
        }
    ]
    
    # Mock read_csv to return DataFrames
    mock_read_csv.return_value = pd.DataFrame({
        "GENE": ["gene1", "gene2"],
        "normZ": [1.5, -1.2],
        "FDR": [0.01, 0.05]
    })
    
    # Mock calc_stats to return statistics
    mock_calc_stats.return_value = {
        "normz_R2_DZ": 0.95,
        "normz_slope_DZ": 0.92,
        "fdr_synth_R2_DZ": 0.88,
        "fdr_synth_slope_DZ": 0.85
    }
    
    # Mock ensure_output_dir to return a simple string path
    mock_ensure_dir.return_value = "/output/dir"
    
    # Call function with plots disabled
    with patch('analysis_pipeline.benchmarking.benchmark_pipeline.Path') as mock_path:
        # Make Path return itself for this test
        mock_path.side_effect = lambda p: p
        
        result = run_direct_comparison(
            file_pairs=file_pairs,
            output_dir="/output/dir",
            generate_plots=False
        )
    
    # Check results
    assert len(result) == 1
    assert result["Experiment"][0] == "experiment1"
    assert result["Contrast"][0] == "contrast1"
    assert result["Analysis"][0] == "DrugZ"
    assert result["normz_R2"][0] == 0.95
    assert result["normz_slope"][0] == 0.92
    
    # Verify calls were made correctly
    mock_ensure_dir.assert_called()
    assert mock_calc_stats.call_count > 0
    mock_to_csv.assert_called_once()


if __name__ == "__main__":
    pytest.main(["-v", __file__]) 