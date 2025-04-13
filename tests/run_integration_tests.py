#!/usr/bin/env python3
"""
Script to run integration tests with the correct Python path.
"""

import os
import sys
import subprocess
from pathlib import Path

# Add the project root to the Python path
# Go up one level from the script's directory (tests/) to get the project root
project_root = Path(__file__).parent.parent.absolute()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
    print(f"Added project root to sys.path: {project_root}")
else:
    print(f"Project root already in sys.path: {project_root}")

# Now run the integration tests
print(f"Running integration tests script: tests/integration/run_all_tests.py")
try:
    result = subprocess.run(
        [sys.executable, "tests/integration/run_all_tests.py"],
        check=False,
        cwd=project_root # Run tests from the project root directory
    )
    sys.exit(result.returncode)
except Exception as e:
    print(f"Error running tests: {e}")
    sys.exit(1) 