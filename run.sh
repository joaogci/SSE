#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-h] [-n <n_threads>] [-o <output_name>] [-t]

-h                  : Help. Shows this text.
-n <n_threads>      : Set number of threads for openMP. Default value: "4"
-o <output_name>    : Relative path + name of output file for simulation. csv file is more convinient. Default value: "output.csv"
-t                  : Test mode. Uses a pre-fixed seed (2)

EOF
    exit 1
}

output_name="output.csv"
n_threads="4"
test=""
cond=""

while getopts chn:o:t opts; do
    case $opts in 
        (c) cond="cond";;
        (h) usage;;
        (n) n_threads=$OPTARG;;
        (o) output_name=$OPTARG;;
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
cd src
make $cond $test
cd ..

echo "[4] - Running the simulation."
echo
./src/main $n_threads tmp/$vtx_name
echo 
