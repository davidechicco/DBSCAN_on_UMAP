#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
set -o xtrace

random_num=$((1 + $RANDOM % 10000))
fileName="../results/output_"$random_num"execution.txt";

Rscript umap_dbscan_EHRs_data.r > $fileName 2> $fileName; 
echo "The end"
