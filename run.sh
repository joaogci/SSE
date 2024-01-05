#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [-c <option>] [-h] [-n <n_threads>] [-m <max_samp>] [-t]

-h                  : Help. Shows this text.
-k <option>         : Compute kinetic coefficients. Options: "full", "ss", "hh", "diag", "offd", "sssh", "hhsh"
-n <n_threads>      : Set number of threads for openMP. Default value: "4"
-m <max_samp>       : Number of maximum samples for conductance
-t                  : Test mode. Uses a pre-fixed seed (2)

EOF
    exit 1
}

n_threads="4"
test=""
kinetic=""
max_samp="1"

while getopts k:hn:m:t opts; do
    case $opts in 
        (k) kinetic=$OPTARG;;
        (h) usage;;
        (n) n_threads=$OPTARG;;
        (m) max_samp=$OPTARG;;
        (t) test="test";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1;;
        (*) usage
    esac
done

if [[ -z "${SSE_DIR}" ]]; then
  echo "The enviroment variable SSE_DIR does not exist."
  echo "Please run the build.sh script on the main directory: "
  echo "source build.sh"
  exit 1 
fi

if [[ ! -f "parameters" ]]; then
  echo "The parameters files does not exist in the current directory. "
  echo "Please copy the file from the main folder, ${SSE_DIR}."
fi

echo "[1] - Generating vertex information."
python3 $SSE_DIR/src/hamiltonian/gen_vtx.py parameters

echo "[2] - Running the simulation."
$SSE_DIR/src/main $n_threads
echo "[2] - Finished the simulation."
