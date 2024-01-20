#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME 

EOF
    exit 1
}

echo "[1] - Exporting directories"
export SSE_DIR=$(pwd)

echo "[2] - Compiling SSE code"
cd src
make
cd ..
echo "[3] -  Successful compilation"

