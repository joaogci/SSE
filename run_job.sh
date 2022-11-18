#!/bin/bash

# Job name:
#SBATCH --job-name=XY_L128
#
# Project:
#SBATCH --account=ec12
#
# Wall time limit:
#SBATCH --time=1-00:00:00
#
# Other parameters:
#SBATCH --mem-per-cpu=3G
#SBATCH --ntasks=1

PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-c] [-h] [-n <n_threads>] [-t]

-c                  : Run conductivity calculations.
-h                  : Help. Shows this text.
-n <n_threads>      : Set number of threads for openMP. Default value: "4"
-t                  : Test mode. Uses a pre-fixed seed (2)

EOF
    exit 1
}

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
set -o pipefail

module --quiet purge  # Reset the modules to the system default
module load GCC/11.2.0
module load GCCcore/11.2.0
module load Python/3.9.6-GCCcore-11.2.0
module load SciPy-bundle/2021.10-foss-2021b
module list

## Do some work:

n_threads="4"
test=""
cond=""

while getopts chn:t opts; do
    case $opts in 
        (c) cond="cond";;
        (h) usage;;
        (n) n_threads=$OPTARG;;
        (t) test="test";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1;;
        (*) usage
    esac
done

echo "[1] - Checking if tmp directory exists."
if [ -d "tmp" ]; then
    echo "tmp directoty exists."
else 
    echo "tmp directory does not exist. Creating directory."
    mkdir tmp
fi

echo "[2] - Generating vertex information."
vtx_name=$(python3 src/vtx/gen_vtx.py read.in)
echo "Saved to file $vtx_name"

echo "[3] - Compiling program."
make $cond $test

echo "[4] - Running the simulation."
echo
./main $n_threads tmp/$vtx_name
echo "[4] - Finished the simulation."
