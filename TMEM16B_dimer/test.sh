#!/bin/bash

for file in ./*_bilayer.inp; do
    base=$(basename "$file")          # e.g. DOPE_bilayer.pdb
    lipid=${base%%_*}                 # e.g. DOPE

    # Check if corresponding .gro exists
    if [[ -f "./${lipid}_bilayer_md3.gro" ]]; then
        continue
    fi

    echo "$lipid"
done
