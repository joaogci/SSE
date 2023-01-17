#!/bin/bash

# Job name:
# no SBATCH --job-name=xy_L128
#
# Project:
#SBATCH --account=ec12
#
# Wall time limit:
# no SBATCH --time=1-00:00:00
#
# Other parameters:
#SBATCH --ntasks=1
# no SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=500M

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
set -o pipefail

module --quiet purge  # Reset the modules to the system default
module load GCC/11.2.0
module load GCCcore/11.2.0
module load Python/3.9.6-GCCcore-11.2.0
module load SciPy-bundle/2021.10-foss-2021b
module load Arb/2.22.1-GCC-11.2.0
module list

./main $1 $2
