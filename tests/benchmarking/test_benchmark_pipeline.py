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
    
    mock_gen_exp2 = MagicMock()
    mock_gen_exp2.name = "experiment2"
    mock_gen_exp2.is_dir.return_value = True
    
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
    
    # Setup published files
    mock_screen1_pdz = MagicMock()
    mock_screen1_pdz.name = "screen1_pDZ.csv"
    mock_screen1_pdz.exists.return_value = True
    
    mock_screen2_pdz = MagicMock()
    mock_screen2_pdz.name = "screen2_pDZ.csv"
    mock_screen2_pdz.exists.return_value = True
    
    mock_screen1_pmgk = MagicMock()
    mock_screen1_pmgk.name = "screen1_pMGK.csv"
    mock_screen1_pmgk.exists.return_value = True
    
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
        elif arg.name.endswith("_pDZ.csv"):
            if arg.name == "screen1_pDZ.csv":
                return mock_screen1_pdz
            elif arg.name == "screen2_pDZ.csv":
                return mock_screen2_pdz
        elif arg.name.endswith("_pMGK.csv"):
            if arg.name == "screen1_pMGK.csv":
                return mock_screen1_pmgk
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
    
    # Call function with plots enabled (default behavior)
    result = run_direct_comparison(
        file_pairs=file_pairs,
        output_dir="/output/dir"
    )
    
    # Check results
    assert len(result) == 1
    assert result["Experiment"][0] == "experiment1"
    assert result["Contrast"][0] == "contrast1"
    assert result["Analysis"][0] == "DrugZ"
    assert result["normz_R2"][0] == 0.95
    assert result["normz_slope"][0] == 0.92
    
    # Verify calls
    mock_ensure_dir.assert_called()
    mock_calc_stats.assert_called_with(
        p_df=mock_read_csv.return_value,
        g_df=mock_read_csv.return_value,
        analysis_type="DZ",
        output_dir=mock_ensure_dir.return_value,
        graph=True  # Plots should be enabled by default
    )
    mock_to_csv.assert_called_once()
    
    # Reset mocks
    mock_ensure_dir.reset_mock()
    mock_calc_stats.reset_mock()
    mock_to_csv.reset_mock()
    
    # Test with plots disabled
    result = run_direct_comparison(
        file_pairs=file_pairs,
        output_dir="/output/dir",
        generate_plots=False
    )
    
    # Verify plots are disabled
    mock_calc_stats.assert_called_with(
        p_df=mock_read_csv.return_value,
        g_df=mock_read_csv.return_value,
        analysis_type="DZ",
        output_dir=mock_ensure_dir.return_value,
        graph=False
    )


if __name__ == "__main__":
    pytest.main(["-v", __file__]) 