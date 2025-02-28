# Installation Guide

This guide walks you through the installation process for the CRISPR Analysis Pipeline.

## Prerequisites

- Python 3.8 or later
- Docker installed and running (for MAGeCK and DrugZ analysis)
- Git (for development installation)

## User Installation

### 1. Clone the repository:

```bash
git clone https://github.com/spejohn/crispr_analysis.git
cd crispr_analysis
```

### 2. Create a virtual environment (recommended):

```bash
# Create a virtual environment
python -m venv venv

# Activate the virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate
```

### 3. Install the package:

```bash
# Install in development mode
pip install -e .

# Or install with all dependencies
pip install -e ".[dev,test]"
```

### 4. Verify the installation:

```bash
# Verify Docker is configured correctly
python -c "from analysis_pipeline.docker.docker_utils import verify_docker; verify_docker()"

# Run a simple test analysis
python -m analysis_pipeline.pipeline --help
```

## Development Installation

For development, you'll need additional dependencies:

```bash
# Install development and test dependencies
pip install -e ".[dev,test]"

# Install pre-commit hooks
pre-commit install
```

## Docker Configuration

### Docker on Windows

1. Install Docker Desktop for Windows with WSL 2 backend
2. Ensure Docker Desktop is running
3. Configure resource sharing for your data directories

### Docker on macOS/Linux

1. Install Docker using your package manager or from Docker's website
2. Start the Docker daemon
3. Ensure your user has permissions to run Docker commands

## Troubleshooting

### Docker Connectivity Issues

If you see errors connecting to Docker:

1. Verify Docker is running: `docker info`
2. Check user permissions: On Linux, ensure your user is in the docker group
3. Restart Docker service

### Path Conversion Issues

On Windows, if you encounter path conversion errors:

1. Use forward slashes in your paths
2. Avoid spaces and special characters in directory names
3. Use relative paths when possible

### Missing Dependencies

If you see ImportError or ModuleNotFoundError:

1. Ensure your virtual environment is activated
2. Reinstall the package: `pip install -e .`
3. Check for missing dependencies: `pip install -r requirements.txt`

## Updating

To update to the latest version:

```bash
git pull
pip install -e .
```

## Uninstallation

```bash
pip uninstall analysis_pipeline
``` 