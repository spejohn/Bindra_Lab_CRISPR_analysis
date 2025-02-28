#!/usr/bin/env python3
"""
Script to build Docker images for CRISPR analysis pipeline.
"""

import os
import subprocess
import argparse
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def build_docker_images(pull: bool = False, no_cache: bool = False) -> bool:
    """
    Build Docker images using docker-compose.
    
    Args:
        pull: Whether to pull base images first
        no_cache: Whether to build without using cache
    
    Returns:
        bool: True if build was successful, False otherwise
    """
    # Get the directory containing this script
    script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
    
    # Change to the docker directory
    os.chdir(script_dir)
    
    logger.info(f"Building Docker images from {script_dir}")
    
    # Build command
    cmd = ["docker-compose", "build"]
    
    if pull:
        cmd.append("--pull")
    
    if no_cache:
        cmd.append("--no-cache")
    
    # Run the command
    try:
        logger.info(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        logger.info("Docker images built successfully")
        logger.debug(result.stdout)
        
        # Verify the images were built
        images = ["crispr-analysis/mageck:latest", "crispr-analysis/drugz:latest"]
        missing_images = []
        
        for image in images:
            img_check = subprocess.run(
                ["docker", "image", "inspect", image],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
            
            if img_check.returncode != 0:
                missing_images.append(image)
        
        if missing_images:
            logger.error(f"Build completed but the following images are missing: {', '.join(missing_images)}")
            return False
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error building Docker images: {e}")
        logger.error(e.stderr)
        return False


def main():
    """
    Main function to build Docker images.
    """
    parser = argparse.ArgumentParser(description="Build Docker images for CRISPR analysis pipeline")
    parser.add_argument("--pull", action="store_true", help="Pull base images before building")
    parser.add_argument("--no-cache", action="store_true", help="Build without using cache")
    
    args = parser.parse_args()
    
    success = build_docker_images(pull=args.pull, no_cache=args.no_cache)
    
    if success:
        print("Docker images built successfully!")
    else:
        print("Failed to build Docker images. See logs for details.")
        exit(1)


if __name__ == "__main__":
    main() 