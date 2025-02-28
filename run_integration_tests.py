#!/usr/bin/env python3
"""
Script to run integration tests with the correct Python path.
"""

import os
import sys
import subprocess
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent.absolute()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Check if package is installed properly
try:
    import analysis_pipeline
    print(f"analysis_pipeline is installed: version {analysis_pipeline.__version__}")
except ImportError:
    print("analysis_pipeline module is not properly installed.")
    print("You may need to run `python install_dev.py` to install it in development mode.")
    print("Continuing anyway with current Python path...")

# Now run the integration tests
print(f"Running integration tests with Python path: {sys.path[0]}")
try:
    result = subprocess.run(
        [sys.executable, "tests/integration/run_all_tests.py"],
        check=False,
        cwd=project_root
    )
    sys.exit(result.returncode)
except Exception as e:
    print(f"Error running tests: {e}")
    sys.exit(1) 