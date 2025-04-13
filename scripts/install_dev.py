#!/usr/bin/env python3
"""
Script to install the CRISPR Analysis Pipeline in development mode.
This makes it easier to test and develop the code.
"""

import os
import sys
import subprocess
from pathlib import Path

def main():
    """
    Install the package in development mode.
    """
    # Get the current directory
    project_root = Path(__file__).parent.absolute()
    
    print(f"Installing CRISPR Analysis Pipeline in development mode...")
    print(f"Project root: {project_root}")
    
    # Install the package in development mode
    cmd = [sys.executable, "-m", "pip", "install", "-e", "."]
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            cwd=project_root,
            capture_output=True,
            text=True
        )
        print("Output:")
        print(result.stdout)
        
        print("\nInstallation successful!")
        print("\nYou can now import the package with:")
        print("  import analysis_pipeline")
        
        # Verify the installation
        print("\nVerifying installation...")
        verify_cmd = [
            sys.executable, 
            "-c", 
            "import analysis_pipeline; print(f'Successfully imported analysis_pipeline version {analysis_pipeline.__version__}')"
        ]
        verify_result = subprocess.run(
            verify_cmd,
            check=False,
            capture_output=True,
            text=True
        )
        
        if verify_result.returncode == 0:
            print(verify_result.stdout.strip())
            return 0
        else:
            print("Verification failed:")
            print(verify_result.stderr)
            return 1
            
    except subprocess.CalledProcessError as e:
        print("Error during installation:")
        print(e.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main()) 