"""
Utility functions for the CRISPR analysis pipeline.
"""

import os
import sys
import time
import datetime
import platform
import functools
import logging
from typing import Optional, Dict, List, Any, Callable, Union
import subprocess
import psutil

# Retry decorator for operations that might fail temporarily
def retry_operation(max_attempts: int = 3, delay: int = 5, 
                   error_types: tuple = (IOError, OSError, ConnectionError)):
    """
    Decorator that retries a function if it fails with specified error types.
    
    Args:
        max_attempts: Maximum number of retry attempts
        delay: Delay between retries in seconds
        error_types: Tuple of exception types to catch and retry
        
    Returns:
        Decorated function that will retry on failure
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            last_exception = None
            
            while attempts < max_attempts:
                try:
                    return func(*args, **kwargs)
                except error_types as e:
                    attempts += 1
                    last_exception = e
                    
                    if attempts == max_attempts:
                        logging.error(f"Operation failed after {max_attempts} attempts: {e}")
                        raise
                        
                    logging.warning(f"Operation failed (attempt {attempts}/{max_attempts}), "
                                   f"retrying in {delay} seconds: {e}")
                    time.sleep(delay)
            
            # This shouldn't be reached due to the raise above, but just in case
            if last_exception:
                raise last_exception
            
        return wrapper
    return decorator


def check_docker_available() -> bool:
    """
    Check if Docker is available and working.
    
    Returns:
        True if Docker is available, False otherwise
    """
    try:
        # Try to import the docker module
        import docker
        
        # Try to connect to the Docker daemon
        client = docker.from_env()
        client.ping()  # Will raise an exception if Docker is not running
        logging.info("Docker is available")
        return True
    except ImportError:
        logging.warning("Docker Python package is not installed")
        return False
    except Exception as e:
        logging.warning(f"Docker is not available: {e}")
        return False


def run_subprocess(command: List[str], timeout: Optional[int] = None) -> Dict[str, Any]:
    """
    Run a subprocess command and return the results.
    
    Args:
        command: Command to run as a list of strings
        timeout: Timeout in seconds
        
    Returns:
        Dictionary with stdout, stderr, and return code
    """
    try:
        # Start the process
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        # Wait for the process to complete
        stdout, stderr = process.communicate(timeout=timeout)
        
        return {
            'stdout': stdout,
            'stderr': stderr,
            'returncode': process.returncode
        }
    except subprocess.TimeoutExpired:
        process.kill()
        stdout, stderr = process.communicate()
        logging.error(f"Command {' '.join(command)} timed out after {timeout} seconds")
        return {
            'stdout': stdout,
            'stderr': stderr,
            'returncode': -1,
            'error': 'Timeout'
        }
    except Exception as e:
        logging.error(f"Error running command {' '.join(command)}: {e}")
        return {
            'stdout': '',
            'stderr': str(e),
            'returncode': -1,
            'error': str(e)
        }


def get_system_info() -> Dict[str, str]:
    """
    Get information about the system.
    
    Returns:
        Dictionary with system information
    """
    info = {
        'platform': platform.platform(),
        'python_version': sys.version,
        'hostname': platform.node(),
        'date': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }
    
    try:
        # Try to get CPU info
        import multiprocessing
        info['cpu_count'] = str(multiprocessing.cpu_count())
    except:
        info['cpu_count'] = 'Unknown'
    
    try:
        # Try to get memory info
        mem = psutil.virtual_memory()
        info['memory_total'] = f"{mem.total / (1024**3):.1f} GB"
        info['memory_available'] = f"{mem.available / (1024**3):.1f} GB"
    except:
        info['memory_total'] = 'Unknown'
        info['memory_available'] = 'Unknown'
    
    return info

def convert_win_path_to_docker(path_str):
    """
    Convert Windows path to Docker-compatible format.
    
    This function converts Windows-style paths to Docker-compatible paths:
    - Replaces backslashes with forward slashes
    - Converts drive letters (e.g., C:) to lowercase Docker format (/c/)
    
    Args:
        path_str (str): Windows path string (e.g., "C:\\path\\to\\dir")
        
    Returns:
        str: Docker-compatible path (e.g., "/c/path/to/dir")
        
    Examples:
        >>> convert_win_path_to_docker("C:\\Users\\test\\data")
        '/c/Users/test/data'
        
        >>> convert_win_path_to_docker("D:\\CRISPR\\results")
        '/d/CRISPR/results'
        
        >>> convert_win_path_to_docker("/home/user/data")  # Non-Windows path
        '/home/user/data'
    """
    import os
    
    if os.name == 'nt':  # Windows
        if ':' in path_str:
            drive, rest = path_str.split(':', 1)
            return f"/{drive.lower()}{rest.replace('\\', '/')}"
    return path_str.replace('\\', '/') 