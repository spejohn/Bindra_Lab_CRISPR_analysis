def check_docker_available() -> bool:
    """
    Check if Docker is available on the system.
    
    Returns:
        bool: True if Docker is available, False otherwise
    """
    import subprocess
    import logging
    
    try:
        result = subprocess.run(
            ["docker", "--version"], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        
        if result.returncode == 0:
            logging.info(f"Docker is available: {result.stdout.strip()}")
            
            # Also check if our images are available
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
                logging.warning(f"Docker is available but the following images are missing: {', '.join(missing_images)}")
                logging.warning("Run 'docker-compose build' in the analysis_pipeline/docker directory to build the images")
                return False
            
            return True
        else:
            logging.warning("Docker is not available or not in PATH")
            return False
    except Exception as e:
        logging.warning(f"Error checking Docker availability: {e}")
        return False 

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