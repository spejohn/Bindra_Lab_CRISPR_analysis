#!/bin/bash
# slurm-submit.sh - Submit script for Snakemake SLURM profile

# Set defaults for common parameters
DEFAULT_PARTITION="day"
DEFAULT_MEM_MB=16000
DEFAULT_TIME_MIN=120

# Parse command line arguments
PARAMS=""
while [ $# -gt 0 ]; do
    case "$1" in
        --cpus-per-task=*)
            CPUS="${1#*=}"
            ;;
        --mem=*)
            MEM="${1#*=}"
            ;;
        --time=*)
            TIME="${1#*=}"
            ;;
        --job-name=*)
            NAME="${1#*=}"
            ;;
        --output=*)
            OUT="${1#*=}"
            ;;
        --error=*)
            ERR="${1#*=}"
            ;;
        --partition=*)
            PARTITION="${1#*=}"
            ;;
        # Add any other parameters you specifically set in config.yaml
        # Default: pass unknown args to sbatch
        *)
            PARAMS="$PARAMS $1"
            ;;
    esac
    shift
done

# Set defaults for missing values
CPUS=${CPUS:-1}
MEM=${MEM:-$DEFAULT_MEM_MB}
TIME=${TIME:-$DEFAULT_TIME_MIN}
PARTITION=${PARTITION:-$DEFAULT_PARTITION}

# Ensure logs directory exists
mkdir -p "$(dirname "$OUT")" "$(dirname "$ERR")"

# Build the sbatch command
sbatch_cmd="sbatch"
sbatch_cmd+=" --parsable"                       # Output only job ID for better parsing
sbatch_cmd+=" --cpus-per-task=$CPUS"            # Set CPU count
sbatch_cmd+=" --mem=${MEM}"                     # Set memory
sbatch_cmd+=" --time=${TIME}"                   # Set time limit
sbatch_cmd+=" --partition=${PARTITION}"         # Set partition
sbatch_cmd+=" --job-name=${NAME}"               # Set job name
sbatch_cmd+=" --output=${OUT}"                  # Set stdout log file
sbatch_cmd+=" --error=${ERR}"                   # Set stderr log file
sbatch_cmd+=" --mail-user=${DEFAULT_MAIL_USER}" # Set email notification
sbatch_cmd+=" --mail-type=${DEFAULT_MAIL_TYPE}" # Set notification type
sbatch_cmd+=" $PARAMS"                          # Add any additional parameters

# Execute the sbatch command
eval "$sbatch_cmd" 