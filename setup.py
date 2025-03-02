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
    name="analysis_pipeline",
    version="0.1.0",
    description="CRISPR Analysis Pipeline - A modular pipeline for CRISPR screening data analysis",
    author="CRISPR Analysis Team",
    packages=find_packages(include=['analysis_pipeline', 'analysis_pipeline.*']),
    install_requires=requirements,
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            'run_pipeline=analysis_pipeline.pipeline:main',
            'run_snakemake=analysis_pipeline.run_snakemake:main',
            'run_qc_analysis=analysis_pipeline.run_qc_analysis:main',
            'run_multiple_screens=analysis_pipeline.run_multiple_screens:main',
            'convert_input_files=analysis_pipeline.convert_input_files:main',
        ],
    },
) 