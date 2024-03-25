#!/bin/bash

# Job name:
#SBATCH --job-name=Dir_R
# Project:
#SBATCH --account=ec12
# Wall time limit:
#SBATCH --time=1-00:00:00
# Other parameters:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=NCPU_R
#SBATCH --mem-per-cpu=750M
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joao.carvalho-inacio@uni-wuerzburg.de


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
set -o pipefail

module --quiet purge  # Reset the modules to the system default
module load GCC/12.3.0
module load GCCcore/12.3.0
module load Python/3.11.3-GCCcore-12.3.0
module load SciPy-bundle/2023.07-gfbf-2023a
module list

$SSE_DIR/src/main NCPU_R

