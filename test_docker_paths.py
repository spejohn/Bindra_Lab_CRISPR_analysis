#!/usr/bin/env python3
"""
Test Docker path conversion functionality.
This is a simple test script that can be run to verify the Docker path conversion works.
"""

import os
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to import the centralized function
try:
    # Add the parent directory to the sys.path to enable imports
    parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if parent_dir not in sys.path:
        sys.path.append(parent_dir)
    
    from analysis_pipeline.core.utils import convert_win_path_to_docker
    logger.info("Successfully imported convert_win_path_to_docker from core.utils")
except ImportError:
    # Fall back to direct import
    try:
        from core.utils import convert_win_path_to_docker
        logger.info("Successfully imported convert_win_path_to_docker from local core.utils")
    except ImportError:
        # Define the function locally as a last resort
        logger.warning("Could not import from core.utils, defining function locally")
        def convert_win_path_to_docker(path_str):
            """
            Convert Windows path to Docker-compatible format.
            
            Args:
                path_str (str): Windows path string (e.g. C:\\path\\to\\dir)
                
            Returns:
                str: Docker-compatible path (e.g. /c/path/to/dir)
            """
            if os.name == 'nt':  # Windows
                if ':' in path_str:
                    drive, rest = path_str.split(':', 1)
                    rest_converted = rest.replace('\\', '/')
                    return f"/{drive.lower()}{rest_converted}"
            return path_str.replace('\\', '/')

def test_path_conversion():
    """
    Test Windows path conversion for Docker compatibility.
    """
    test_paths = [
        "C:\\Users\\test\\data",
        "D:\\CRISPR\\results",
        "/home/user/data",
        "relative/path",
        "C:\\Users\\test\\path with spaces",
        "E:\\Data\\CRISPR\\screens\\experiment_1\\counts.csv",
        "\\\\server\\share\\data",
        "C:"
    ]
    
    logger.info("Testing Docker path conversion:")
    for path in test_paths:
        converted = convert_win_path_to_docker(path)
        logger.info(f"  Original: {path}")
        logger.info(f"  Converted: {converted}")
        logger.info("")
    
    # Test current working directory
    cwd = os.getcwd()
    logger.info(f"Current working directory: {cwd}")
    logger.info(f"Converted: {convert_win_path_to_docker(cwd)}")
    
    # Test modules that use the conversion function
    try:
        # These imports will fail if the modules aren't in the path
        # We're just testing that these imports would work if the path is set up
        logger.info("Testing imports of modules that use path conversion")
        logger.info("Note: These may fail if not running as a package, which is expected")
        
        try:
            # Check if docker package is installed
            try:
                import docker
                logger.info("Docker package is installed")
            except ImportError:
                logger.warning("Docker package is not installed, some functionality may be limited")
                
            # Try to import modules that use the conversion function
            try:
                from analysis_pipeline.analysis.run_drugz_docker import run_drugz_docker
                logger.info("  Successfully imported run_drugz_docker.run_drugz_docker")
            except ImportError as e:
                logger.warning(f"  Could not import run_drugz_docker: {e}")
                
            try:
                from analysis_pipeline.analysis.mageck_analysis import mageck_test_analysis
                logger.info("  Successfully imported mageck_analysis.mageck_test_analysis")
            except ImportError as e:
                logger.warning(f"  Could not import mageck_analysis: {e}")
        except Exception as e:
            logger.warning(f"Error testing module imports: {e}")
    except Exception as e:
        logger.warning(f"Error in import tests: {e}")

if __name__ == "__main__":
    logger.info(f"Operating system: {os.name}")
    logger.info(f"Python version: {sys.version}")
    test_path_conversion()
    logger.info("Path conversion test completed.") 