"""
Provides function to execute commands within containers (Apptainer/Docker).
"""

import subprocess
import logging
import os
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple

# Attempt to import Docker utilities, handle potential import errors gracefully
try:
    # Assumes core is accessible from the execution context
    from docker.docker_utils import run_docker_container
    DOCKER_UTILS_AVAILABLE = True
except ImportError:
    DOCKER_UTILS_AVAILABLE = False
    # Define a placeholder if import fails, so the function signature exists
    def run_docker_container(*args, **kwargs) -> Tuple[int, str]:
        logger = logging.getLogger(__name__)
        logger.error("Docker utilities not available. Cannot run via Docker.")
        return (1, "Docker utilities not available.")

logger = logging.getLogger(__name__)


def execute_in_container(
    command_list: List[str],
    container_image: str,
    mount_map: Dict[str, str],
    working_dir: Optional[str] = None,
    use_apptainer: bool = True,
    timeout: int = 3600,
    stream_logs: bool = True
) -> Tuple[int, str]:
    """
    Executes a command list within a specified container using Apptainer or Docker.

    Prioritizes Apptainer if use_apptainer is True.

    Args:
        command_list: The command and its arguments as a list.
        container_image: Path to the Apptainer SIF file or Docker image URI (e.g., docker://ubuntu:latest).
        mount_map: Dictionary mapping host paths to container paths for mounting.
        working_dir: Optional path inside the container to set as the working directory.
        use_apptainer: If True, attempt to use Apptainer first. If False or Apptainer fails,
                       fallback to Docker (if available).
        timeout: Maximum execution time in seconds.
        stream_logs: If True, print captured output after successful execution.

    Returns:
        Tuple containing the exit code (int) and captured output (str: stdout + stderr).
    """

    if use_apptainer:
        logger.info(f"Attempting execution via Apptainer for image {container_image}")
        apptainer_cmd = ["apptainer", "exec"]
        bind_args = []
        for host_path, container_path in mount_map.items():
            abs_host_path = os.path.abspath(host_path)
            # Create host directory if it doesn't exist to avoid Apptainer errors
            Path(abs_host_path).mkdir(parents=True, exist_ok=True)
            bind_args.extend(["--bind", f"{abs_host_path}:{container_path}"])

        apptainer_cmd.extend(bind_args)

        if working_dir:
            # Ensure working_dir is an absolute path within the container context
            if not os.path.isabs(working_dir):
                 logger.warning(f"Provided working_dir '{working_dir}' is not absolute. Apptainer might interpret it relative to its default.")
            apptainer_cmd.extend(["--pwd", working_dir])

        apptainer_cmd.append(container_image)
        apptainer_cmd.extend(command_list)

        logger.info(f"Running Apptainer command: {' '.join(apptainer_cmd)}")
        try:
            result = subprocess.run(
                apptainer_cmd,
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout
            )
            output = f"APPTAINER STDOUT:\n{result.stdout}\nAPPTAINER STDERR:\n{result.stderr}"
            if result.returncode == 0:
                logger.info(f"Apptainer command successful (Exit Code: {result.returncode}).")
                if stream_logs:
                    print(output)
                return result.returncode, output
            else:
                logger.error(f"Apptainer command failed (Exit Code: {result.returncode}) for image {container_image}")
                logger.error(output)
                # Decide whether to fallback to Docker here or let the caller handle it
                # For now, return the Apptainer failure
                return result.returncode, output

        except FileNotFoundError:
            logger.error("Apptainer command not found. Is Apptainer installed and in PATH?")
            # Fallback or error? For now, error out if Apptainer was specifically requested/prioritized
            return 1, "Apptainer command not found."
        except subprocess.TimeoutExpired:
            logger.error(f"Apptainer command timed out after {timeout} seconds.")
            return 1, f"Command timed out after {timeout} seconds."
        except Exception as e:
            logger.error(f"Error running Apptainer command: {e}")
            return 1, f"Error running Apptainer command: {e}"

    # If not using Apptainer or if it failed (decision needed above) - try Docker
    logger.info(f"Attempting execution via Docker for image {container_image}")
    if not DOCKER_UTILS_AVAILABLE:
        logger.error("Docker execution requested, but Docker utilities are not available.")
        return 1, "Docker utilities not available."

    # Convert mount_map to Docker SDK volumes format
    volumes_dict = {}
    for host_path, container_path in mount_map.items():
        abs_host_path = os.path.abspath(host_path)
        Path(abs_host_path).mkdir(parents=True, exist_ok=True)
        volumes_dict[abs_host_path] = {'bind': container_path, 'mode': 'rw'}

    # *** TODO: Modify docker_utils.run_docker_container to accept working_dir ***
    if working_dir:
        logger.warning("Docker execution path needs 'run_docker_container' to be adapted for 'working_dir'. Ignoring working_dir for Docker execution currently.")

    # Convert Docker image URI if needed (e.g., docker://ubuntu -> ubuntu)
    docker_image_name = container_image
    if docker_image_name.startswith("docker://"):
        docker_image_name = docker_image_name[len("docker://"):]

    return run_docker_container(
        command=command_list,
        volumes=volumes_dict,
        image=docker_image_name, # Use potentially modified name
        timeout=timeout,
        stream_logs=stream_logs,
        working_dir=working_dir
    ) 