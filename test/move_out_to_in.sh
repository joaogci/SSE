#!/bin/bash
shopt -s extglob
 
if [ -e confout_0 ]; then
  for file_out in confout_+([0-9]); do
    file_in=${file_out/out/in}
    echo "mv $file_out to $file_in"
    mv "$file_out" "$file_in"
  done
fi
