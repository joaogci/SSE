#!/bin/bash
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
