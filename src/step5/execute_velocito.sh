#!/bin/bash

# This script aims to execute the velocyto program on LysoDC data.
# It has to be used in a Docker container issued from jpglab_lysodc_rnavelocity docker image
# 

# Provide the folder where the 10x data are stored
export RAW_DATA_PATH=${WORKING_DIR}/data/raw/DIR180504_NB501865_0090_AHM32WBGX5/10635173
# Provide the path to the genome annotation file
export GENOME_ANNOTATION_PATH=${WORKING_DIR}/data/annotation/genome/grch38/Mus_musculus.GRCm38.90.gtf
# Provide the folder where to put the Velocyto results (loom file)
export OUTPUT_PATH=${WORKING_DIR}/data/step5/output
mkdir -p $OUTPUT_PATH

# Ensure Samtools is in the PATH
export PATH=/samtools/bin:$PATH

# Execute Velocyto on 10x data
cd $OUTPUT_PATH
velocyto run10x $RAW_DATA_PATH $GENOME_ANNOTATION_PATH

# Move the Velocyto results (put by the tool in the input folder) in the desired output folder
mv $RAW_DATA_PATH/velocyto/* $OUTPUT_PATH
rmdir $RAW_DATA_PATH/velocyto
