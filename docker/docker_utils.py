"""
Docker utility functions for CRISPR analysis pipeline.
"""

import time
import logging
import docker
from pathlib import Path
from typing import Dict, List, Optional, Union, Any, Tuple

# Import Docker settings from config
from analysis_pipeline.core.config import (
    DOCKER_IMAGE,
    DOCKER_RESTART_POLICY,
    DOCKER_RUN_TIMEOUT
)

# Initialize Docker client
try:
    DOCKER_CLIENT = docker.from_env()
    logging.info("Docker client initialized successfully")
except Exception as e:
    logging.error(f"Error initializing Docker client: {e}")
    DOCKER_CLIENT = None


def verify_docker() -> bool:
    """
    Verify that Docker is available and the MAGeCK image exists.
    
    Returns:
        True if Docker is available and the image exists, False otherwise
    """
    if DOCKER_CLIENT is None:
        logging.error("Docker client not initialized")
        return False
    
    try:
        # Check Docker version
        version = DOCKER_CLIENT.version()
        logging.info(f"Docker version: {version.get('Version', 'unknown')}")
        
        # Check if image exists
        images = DOCKER_CLIENT.images.list(name=DOCKER_IMAGE)
        if images:
            logging.info(f"Docker image found: {DOCKER_IMAGE} ({images[0].tags})")
            return True
        else:
            logging.warning(f"Docker image not found: {DOCKER_IMAGE}")
            return pull_docker_image(DOCKER_IMAGE)
            
    except Exception as e:
        logging.error(f"Error verifying Docker: {e}")
        return False


def pull_docker_image(image_name: str) -> bool:
    """
    Pull the specified Docker image.
    
    Args:
        image_name: Name of the Docker image to pull
        
    Returns:
        True if successful, False otherwise
    """
    if DOCKER_CLIENT is None:
        logging.error("Docker client not initialized")
        return False
    
    try:
        logging.info(f"Pulling Docker image: {image_name}")
        DOCKER_CLIENT.images.pull(image_name)
        logging.info(f"Successfully pulled Docker image: {image_name}")
        return True
    except Exception as e:
        logging.error(f"Error pulling Docker image {image_name}: {e}")
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
        logging.error(error_msg)
        return (1, error_msg)
    
    container = None
    
    try:
        # Log the command and volumes
        if isinstance(command, list):
            cmd_str = " ".join(map(str, command))
        else:
            cmd_str = command
            
        logging.info(f"Running Docker command: {cmd_str}")
        for host_path, container_info in volumes.items():
            logging.info(f"  {host_path} -> {container_info['bind']} ({container_info['mode']})")
        
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
                logging.info(f"Container output: {log_line}")
        
        # Wait for the container to finish with timeout
        start_time = time.time()
        while container.status != "exited" and time.time() - start_time < timeout:
            container.reload()  # Update container status
            time.sleep(1)
            
        # Check if the container timed out
        if container.status != "exited":
            logging.error(f"Container timed out after {timeout} seconds")
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
            logging.error(error_msg)
        else:
            logging.info(f"Container completed successfully (exit code: {exit_code})")
        
        # Clean up container
        try:
            container.remove(force=True)
        except Exception as e:
            logging.warning(f"Could not remove container: {e}")
        
        return (exit_code, logs)
        
    except Exception as e:
        error_msg = f"Error running Docker container: {str(e)}"
        logging.error(error_msg)
        
        # Try to get logs and clean up if container was created
        if container:
            try:
                logs = container.logs().decode()
                logging.error(f"Container logs:\n{logs}")
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
        logging.error("Docker client not initialized")
        return False
    
    try:
        # For the MAGeCK image, test the mageck version command
        test_command = ["mageck", "--version"]
        logging.info(f"Testing Docker container with command: {' '.join(test_command)}")
        
        exit_code, output = run_docker_container(
            command=test_command,
            volumes={},  # No volumes needed for this test
            image=image,
            timeout=30,  # Short timeout for a simple command
            stream_logs=False
        )
        
        if exit_code == 0:
            logging.info(f"Docker container test successful: {output.strip()}")
            return True
        else:
            logging.error(f"Docker container test failed: {output}")
            return False
            
    except Exception as e:
        logging.error(f"Error testing Docker container: {e}")
        return False 