#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-b] [-i <input_name>] [-h] [-n <n_threads>] [-o <output_name>] [-s <series_name>] [-t]

-b                  : Use Heat Bath for transition probabilities (slower method)
-i <input_name>     : Relative path + name of input file for simulation. Default value: "input.txt"
-h                  : Help. Shows this text.
-n <n_threads>      : Set number of threads for openMP. Default value: "4"
-o <output_name>    : Relative path + name of output file for simulation. csv file is more convinient. Default value: "output.csv"
-s <series_name>    : Relative path + name of correlation series file for simluation. csv is more convinient.
-t                  : Test mode. Uses a pre-fixed seed (2) and runs with the test temperatures (1/T = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0}).

EOF
    exit 1
}

input_name="input.txt"
output_name="output.csv"
series_name=""
n_threads="4"
test=""
hb=0

while getopts bi:hn:o:s:t opts; do
    case $opts in 
        (b) hb=1;;
        (i) input_name=$OPTARG;;
        (h) usage;;
        (n) n_threads=$OPTARG;;
        (o) output_name=$OPTARG;;
        (s) series_name=$OPTARG; series="series"; n_threads="1";;
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
vtx_name=$(python3 src/vtx/gen_vtx.py $input_name $hb)
echo "Saved to file $vtx_name"

echo "[3] - Compiling program."
cd src
if [[ "$series" == "series" ]]; then
    make $series
else
    make $test
fi

cd ..
echo "[4] - Running the simulation."
echo
./src/main $n_threads $input_name tmp/$vtx_name $output_name $series_name 
echo 

echo "[5] - Removing vertex file."
rm tmp/${vtx_name}

