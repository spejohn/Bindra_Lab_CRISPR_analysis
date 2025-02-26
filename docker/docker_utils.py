"""
Docker utility functions for CRISPR analysis pipeline.
"""

import time
import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union, Any, Tuple

# Configure logging if not already configured
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to import Docker
DOCKER_AVAILABLE = False
DOCKER_CLIENT = None
try:
    import docker
    DOCKER_AVAILABLE = True
    logger.info("Docker Python package is installed")
except ImportError:
    logger.warning("Docker Python package not installed. Docker functionality will be limited.")
    logger.warning("You may need to run: pip install docker")

# Import Docker settings from config
try:
    # Try relative import for package usage
    from ..core.config import (
        DOCKER_IMAGE,
        DOCKER_RESTART_POLICY,
        DOCKER_RUN_TIMEOUT
    )
    from ..core.utils import convert_win_path_to_docker
except ImportError:
    # Fall back to direct import for direct script execution
    try:
        from core.config import (
            DOCKER_IMAGE,
            DOCKER_RESTART_POLICY,
            DOCKER_RUN_TIMEOUT
        )
        from core.utils import convert_win_path_to_docker
    except ImportError:
        # Default values if config cannot be imported
        DOCKER_IMAGE = "spejohn/mageck:latest"
        DOCKER_RESTART_POLICY = {"Name": "no"}
        DOCKER_RUN_TIMEOUT = 3600
        
        # Define the conversion function directly if we can't import it
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
                    return f"/{drive.lower()}{rest.replace('\\', '/')}"
            return path_str.replace('\\', '/')

# Initialize Docker client using multiple fallback approaches
if DOCKER_AVAILABLE:
    try:
        # Try multiple approaches to initialize the Docker client
        # based on different versions of the docker package
        docker_initialized = False
        
        # Approach 1: Try from_env() directly
        try:
            DOCKER_CLIENT = docker.from_env()
            logger.info("Docker client initialized successfully using docker.from_env()")
            docker_initialized = True
        except (AttributeError, TypeError) as e:
            logger.warning(f"Could not initialize Docker client using docker.from_env(): {e}")
        
        # Approach 2: Try DockerClient.from_env()
        if not docker_initialized:
            try:
                if hasattr(docker, 'DockerClient'):
                    DOCKER_CLIENT = docker.DockerClient.from_env()
                    logger.info("Docker client initialized successfully using DockerClient.from_env()")
                    docker_initialized = True
            except (AttributeError, TypeError) as e:
                logger.warning(f"Could not initialize Docker client using DockerClient.from_env(): {e}")
        
        # Approach 3: Try client.from_env()
        if not docker_initialized:
            try:
                if hasattr(docker, 'client') and hasattr(docker.client, 'from_env'):
                    DOCKER_CLIENT = docker.client.from_env()
                    logger.info("Docker client initialized successfully using client.from_env()")
                    docker_initialized = True
            except (AttributeError, TypeError) as e:
                logger.warning(f"Could not initialize Docker client using client.from_env(): {e}")
        
        # Approach 4: Try APIClient
        if not docker_initialized:
            try:
                if hasattr(docker, 'APIClient'):
                    DOCKER_CLIENT = docker.APIClient()
                    logger.info("Docker client initialized successfully using APIClient")
                    docker_initialized = True
            except (AttributeError, TypeError) as e:
                logger.warning(f"Could not initialize Docker client using APIClient: {e}")
        
        if not docker_initialized:
            logger.error("All Docker client initialization methods failed")
    except Exception as e:
        logger.error(f"Unhandled error initializing Docker client: {e}")

# Define function to check if Docker CLI is available
def check_docker_cli_available() -> bool:
    """
    Check if Docker CLI is available on the system.
    
    Returns:
        bool: True if Docker CLI is available, False otherwise
    """
    try:
        result = subprocess.run(
            ["docker", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        if result.returncode == 0:
            logger.info(f"Docker CLI available: {result.stdout.strip()}")
            return True
        else:
            logger.warning("Docker CLI not available or not in PATH")
            return False
    except Exception as e:
        logger.warning(f"Error checking Docker CLI availability: {e}")
        return False


def verify_docker() -> bool:
    """
    Verify that Docker is available and the MAGeCK image exists.
    
    Returns:
        True if Docker is available and the image exists, False otherwise
    """
    # Check if Docker CLI is available (as a fallback)
    docker_cli_available = check_docker_cli_available()
    
    if not DOCKER_AVAILABLE:
        logger.error("Docker Python package not installed")
        if docker_cli_available:
            logger.info("Docker CLI is available, but Python package is missing")
        return False
    
    if DOCKER_CLIENT is None:
        logger.error("Docker client not initialized")
        if docker_cli_available:
            logger.info("Docker CLI is available, but Python client initialization failed")
        return False
    
    try:
        # Try to check Docker version
        try:
            version = DOCKER_CLIENT.version()
            logger.info(f"Docker version: {version.get('Version', 'unknown')}")
        except Exception as e:
            logger.warning(f"Could not get Docker version: {e}")
            if docker_cli_available:
                logger.info("Using Docker CLI to check version instead")
                return docker_cli_available
            return False
        
        # Try to check if image exists
        try:
            images = DOCKER_CLIENT.images.list(name=DOCKER_IMAGE)
            if images:
                logger.info(f"Docker image found: {DOCKER_IMAGE}")
                return True
            else:
                logger.warning(f"Docker image not found: {DOCKER_IMAGE}")
                return pull_docker_image(DOCKER_IMAGE)
        except Exception as e:
            logger.warning(f"Could not check for Docker image: {e}")
            if docker_cli_available:
                logger.info("Using Docker CLI to check for image instead")
                return docker_cli_available
            return False
            
    except Exception as e:
        logger.error(f"Error verifying Docker: {e}")
        return docker_cli_available  # Fall back to CLI if available


def pull_docker_image(image_name: str) -> bool:
    """
    Pull the specified Docker image.
    
    Args:
        image_name: Name of the Docker image to pull
        
    Returns:
        True if successful, False otherwise
    """
    if DOCKER_CLIENT is None:
        logger.error("Docker client not initialized")
        return False
    
    try:
        logger.info(f"Pulling Docker image: {image_name}")
        DOCKER_CLIENT.images.pull(image_name)
        logger.info(f"Successfully pulled Docker image: {image_name}")
        return True
    except Exception as e:
        logger.error(f"Error pulling Docker image {image_name}: {e}")
        return False


def run_docker_container(
    command: Union[str, List[str]],
    volumes: Dict[str, Dict[str, str]],
    image: str = DOCKER_IMAGE,
    timeout: int = DOCKER_RUN_TIMEOUT,
    stream_logs: bool = True
) -> Tuple[int, str]:
    """
    Run a Docker container with the specified command and volumes.
    
    Args:
        command: Command to run in the container
        volumes: Dictionary mapping host paths to container paths
        image: Docker image to use
        timeout: Timeout in seconds
        stream_logs: Whether to stream logs to the console
        
    Returns:
        Tuple of (exit_code, output)
    """
    if DOCKER_CLIENT is None:
        error_msg = "Docker client not initialized"
        logger.error(error_msg)
        return (1, error_msg)
    
    container = None
    
    try:
        # Log the command and volumes
        if isinstance(command, list):
            cmd_str = " ".join(map(str, command))
        else:
            cmd_str = command
            
        logger.info(f"Running Docker command: {cmd_str}")
        
        # Convert Windows paths for Docker if needed
        if os.name == 'nt':
            windows_volumes = {}
            for host_path, container_info in volumes.items():
                docker_path = convert_win_path_to_docker(host_path)
                windows_volumes[docker_path] = container_info
                logger.info(f"  Windows path conversion: {host_path} -> {docker_path}")
            volumes = windows_volumes
        
        # Log the final volume mappings
        for host_path, container_info in volumes.items():
            logger.info(f"  {host_path} -> {container_info['bind']} ({container_info['mode']})")
        
        # Create and start the container
        container = DOCKER_CLIENT.containers.run(
            image=image,
            command=command,
            volumes=volumes,
            detach=True,
            remove=False,  # Keep container for inspection on error
            restart_policy=DOCKER_RESTART_POLICY
        )
        
        # Stream logs if requested
        if stream_logs:
            log_stream = container.logs(stream=True, follow=True)
            for line in log_stream:
                log_line = line.decode().strip()
                logger.info(f"Container output: {log_line}")
        
        # Wait for the container to finish with timeout
        start_time = time.time()
        while container.status != "exited" and time.time() - start_time < timeout:
            container.reload()  # Update container status
            time.sleep(1)
            
        # Check if the container timed out
        if container.status != "exited":
            logger.error(f"Container timed out after {timeout} seconds")
            container.stop()
            logs = container.logs().decode()
            try:
                container.remove(force=True)
            except:
                pass
            return (1, f"Container timed out. Logs:\n{logs}")
        
        # Get container logs
        logs = container.logs().decode()
        
        # Check exit code
        result = container.wait()
        exit_code = result["StatusCode"]
        
        if exit_code != 0:
            error_msg = f"Container exited with code {exit_code}. Error: {result.get('Error', '')}\nLogs:\n{logs}"
            logger.error(error_msg)
        else:
            logger.info(f"Container completed successfully (exit code: {exit_code})")
        
        # Clean up container
        try:
            container.remove(force=True)
        except Exception as e:
            logger.warning(f"Could not remove container: {e}")
        
        return (exit_code, logs)
        
    except Exception as e:
        error_msg = f"Error running Docker container: {str(e)}"
        logger.error(error_msg)
        
        # Try to get logs and clean up if container was created
        if container:
            try:
                logs = container.logs().decode()
                logger.error(f"Container logs:\n{logs}")
                container.remove(force=True)
            except:
                logs = "Could not retrieve logs"
        else:
            logs = "Container was not created"
            
        return (1, f"{error_msg}\n{logs}")


def test_docker_container(image: str = DOCKER_IMAGE) -> bool:
    """
    Test a simple command to verify the Docker container works.
    
    Args:
        image: Docker image to test
        
    Returns:
        True if successful, False otherwise
    """
    if DOCKER_CLIENT is None:
        logger.error("Docker client not initialized")
        return False
    
    try:
        # For the MAGeCK image, test the mageck version command
        test_command = ["mageck", "--version"]
        logger.info(f"Testing Docker container with command: {' '.join(test_command)}")
        
        exit_code, output = run_docker_container(
            command=test_command,
            volumes={},  # No volumes needed for this test
            image=image,
            timeout=30,  # Short timeout for a simple command
            stream_logs=False
        )
        
        if exit_code == 0:
            logger.info(f"Docker container test successful: {output.strip()}")
            return True
        else:
            logger.error(f"Docker container test failed: {output}")
            return False
            
    except Exception as e:
        logger.error(f"Error testing Docker container: {e}")
        return False 