#!/bin/bash
# Usage: ./juno_nextflow.sh run pipelines/main.nf
module load nextflow
module load java/jdk-11.0.11
module load singularity/3.7.1

SCRIPT_DIR=$(dirname "$0")
echo $SCRIPT_DIR

# input_dir = "/work/greenbaum/users/suns3/MSKCC/iacobuzio/organoids_PDAC/ALN/REDITOOLS_OUTPUT/SINE"
# output_dir = "/work/greenbaum/users/lih7/test/REDItools2_test/analysis/rna-editing-analysis/outputs"

DEFAULT_CONFIG=$SCRIPT_DIR"/nextflow.config"
read -e -p "Config file path [${DEFAULT_CONFIG}]: " CONFIG_PATH

if [ -z "${CONFIG_PATH}" ]; then
  CONFIG_PATH="${DEFAULT_CONFIG}"
fi

read -e -p "Input directory : " input_dir

if [ -z "${input_dir}" ]; then
  echo "invalid input path!"
fi

DEFAULT_OUTPUT_DIR='results'
read -e -p "Output directory [$DEFAULT_OUTPUT_DIR]: " output_dir

if [ -z "${output_dir}" ]; then
  output_dir=$DEFAULT_OUTPUT_DIR
fi
mkdir -p $output_dir
output_dir=$(realpath $output_dir)

DEFAULT_EMAIL=$USER'@mskcc.org'
read -e -p "Notify this Email Address [$DEFAULT_EMAIL]: " EMAIL_ADDRESS
if [ -z "${EMAIL_ADDRESS}" ]; then
  EMAIL_ADDRESS=$DEFAULT_EMAIL
fi

printf "Running RNAediting analysis pipeline:\n----------------\nInput directory: %s\nOutput directory: %s\n----------------\n" "$input_dir" "$output_dir"

# Set project specfic env variables
mkdir -p ${HOME}/.singularity/cache/nextflow
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:=${HOME}/.singularity/cache/nextflow}"

# Run nextflow and pass through args
nohup nextflow run $SCRIPT_DIR/dags/rna_editing.nf -c "${CONFIG_PATH}" -w $SCRIPT_DIR/work --launchDir $SCRIPT_DIR --input_dir $input_dir --output_dir $output_dir -N $EMAIL_ADDRESS -profile juno > run.log 2>&1 &
