#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-i <input_name>] [-h] [-n <n_threads>] [-o <output_name>] [-t]

-b                  : Use Heat Bath for transition probabilities (slower method)
-c                  : Clear vtx tmp directory.
-i <input_name>     : Relative path + name of input file for simulation.
-h                  : Help. Shows this text.
-n <n_threads>      : Set number of threads for openMP.
-o <output_name>    : Relative path + name of output file for simulation. csv file is more convinient.
-t                  : Test mode.

EOF
    exit 1
}

input_name="input.txt"
output_name="output.csv"
n_threads="4"
test=""
clear=0
hb=0

while getopts bci:hn:o:t opts; do
    case $opts in 
        (b) hb=1;;
        (c) clear=1;;
        (i) input_name=$OPTARG;;
        (h) usage;;
        (n) n_threads=$OPTARG;;
        (t) test="test";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1;;
        (*) usage
    esac
done

if [ $clear -eq 1 ]; then 
    echo " ------ "
    echo "Cleaning tmp directory."
    rm -r src/vtx/tmp
fi

echo " ------ "
echo "Checking if tmp directory exists."
if [ -d "src/vtx/tmp" ]; then
    echo "tmp directoty exists."
else 
    echo "tmp directory does not exist. Creating directory."
    mkdir src/vtx/tmp
fi
echo " ------ "

echo "Generating vertex information."
vtx_name=$(python3 src/vtx/gen_vtx.py $input_name $hb)
echo "Saved to file $vtx_name"
echo " ------ "

echo "Compiling program."
cd src
make $test
echo " ------ "

echo "Running the simulation."
echo
./main $n_threads ../$input_name vtx/tmp/$vtx_name ../$output_name 



