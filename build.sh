#!/bin/bash
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [option]

option              : local, fox

EOF
}

if [ $# -lt 1 ]; then
  usage
  return 1
fi
if [ ! "$1" = "local" ] && [ ! "$1" = "fox" ]; then
  usage
  return 1 
fi

if [ "$1" = "fox" ]; then
  module --quiet purge  # Reset the modules to the system default
  module load GCC/12.3.0
  module load GCCcore/12.3.0
  module load Python/3.11.3-GCCcore-12.3.0
  module load SciPy-bundle/2023.07-gfbf-2023a
  module list
fi

echo "[1] - Exporting directories"
export SSE_DIR=$(pwd)

echo "[2] - Compiling SSE code"
cd src
make
cd ..
echo "[3] -  Successful compilation"

