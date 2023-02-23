#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-c <option>] [-h] [-n <n_threads>] [-m <max_samp>] [-t]

-c <option>         : Run conductivity calculations. Options: "both", "heat", "spin"
-h                  : Help. Shows this text.
-n <n_threads>      : Set number of threads for openMP. Default value: "4"
-m <max_samp>       : Number of maximum samples for conductance
-t                  : Test mode. Uses a pre-fixed seed (2)

EOF
    exit 1
}

n_threads="4"
test=""
cond=""
max_samp="1"

while getopts c:hn:m:t opts; do
    case $opts in 
        (c) cond=$OPTARG;;
        (h) usage;;
        (n) n_threads=$OPTARG;;
        (m) max_samp=$OPTARG;;
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
./main $n_threads tmp/$vtx_name $max_samp
echo "[4] - Finished the simulation."
