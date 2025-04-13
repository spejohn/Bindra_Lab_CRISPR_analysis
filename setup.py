#!/usr/bin/env python3
"""
Setup script for CRISPR Analysis Pipeline.
This script configures the package for installation with pip.
"""

from setuptools import setup, find_packages
import os
import sys

# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Read requirements from requirements.txt
with open(os.path.join(current_dir, 'requirements.txt')) as f:
    requirements = f.read().splitlines()

setup(
    name="crispr_analysis_pipeline",
    version="0.2.0",
    description="CRISPR Analysis Pipeline - A modular pipeline for CRISPR screening data analysis",
    author="Spenser Johnson",
    # Find packages automatically, excluding tests, scripts, docs, etc.
    packages=find_packages(exclude=['tests*', 'scripts*', 'docs*', 'test_output*', '*egg-info*']),
    install_requires=requirements,
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            # Remove analysis_pipeline prefix and the non-existent run_pipeline
            'run_snakemake=run_snakemake:main',
            'run_qc_analysis=run_qc_analysis:main',
            'run_multiple_screens=run_multiple_screens:main',
        ],
    },
) 