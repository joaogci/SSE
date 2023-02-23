#!/bin/bash

module use $HOME/.local/easybuild/modules/all
module --quiet purge  # Reset the modules to the system default
module load GCC/11.2.0
module load GCCcore/11.2.0
module load Python/3.9.6-GCCcore-11.2.0
module load SciPy-bundle/2021.10-foss-2021b
module load Arb/2.22.1-GCC-11.2.0

PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-c] [-h] [-n <n_threads>] [-m <max_samp>] [-t]

-c <option>         : Run conductivity calculations. Options: "both", "heat", "spin"
-h                  : Help. Shows this text.
-j <job_name>       : Job name.
-n <n_threads>      : Set number of threads for openMP. Default value: "4"
-m <max_samp>       : Number of maximum samples for conductance
-t                  : Test mode. Uses a pre-fixed seed (2)
-r <DD-HH:MM:SS>    : Time for the job to run.

EOF
    exit 1
}

n_threads="4"
test=""
cond=""
job_name=""
time=""
max_samp="1"

while getopts c:hj:n:m:tr: opts; do
    case $opts in 
        (c) cond=$OPTARG;;
        (h) usage;;
        (j) job_name=$OPTARG;;
        (n) n_threads=$OPTARG;;
        (m) max_samp=$OPTARG;;
        (t) test="test";;
        (r) time=$OPTARG;;
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

echo "[4] - Moving files to another directory."
mkdir $job_name
cp read.in $job_name
cp matsubara.in $job_name
cp beta.in $job_name
cp submit_job.sh $job_name
mv main $job_name
mv tmp/$vtx_name $job_name
cd $job_name

echo "[5] - Submitting job."
echo
sbatch -c $n_threads -J $job_name -t $time ./submit_job.sh $n_threads $vtx_name $max_samp
