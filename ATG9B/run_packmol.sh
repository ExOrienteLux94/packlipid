#!/bin/bash

for file in ./*.inp; do
    fname=$(basename "$file" .inp)
    prefix=${fname%%_*}

    if [[ -f "./${prefix}.pdb" ]]; then
        continue
    fi

    packmol < "$file"

done
