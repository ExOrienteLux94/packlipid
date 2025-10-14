#!/bin/bash

for file in ./*.inp; do
    fname=$(basename "$file" .inp)
    prefix=${fname%%_*}

    if [[ -f "./${prefix}_bilayer_md3.gro" ]]; then
        continue
    fi

    echo "Processing $prefix ..."
    gmx editconf -f "${prefix}.pdb" -o "${prefix}_bilayer.gro"
    echo -e "del 2 - 20\n r SOL \n 0 & ! 1 & ! r SOL \n name 2 SOL \n name 3 MEMB \n q\n" | gmx make_ndx -f "${prefix}_bilayer.gro" -o "${prefix}_bilayer.ndx"
    gmx grompp -f minim0.mdp -c "${prefix}_bilayer.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_minim0.tpr" -maxwarn 1
    gmx mdrun -deffnm "${prefix}_bilayer_minim0"
    gmx grompp -f minim1.mdp -c "${prefix}_bilayer_minim0.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_minim1.tpr" -maxwarn 1
    gmx mdrun -deffnm "${prefix}_bilayer_minim1"
    gmx grompp -f minim2.mdp -c "${prefix}_bilayer_minim1.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_minim2.tpr" -maxwarn 1
    gmx mdrun -deffnm "${prefix}_bilayer_minim2"
    gmx grompp -f md1.mdp -c "${prefix}_bilayer_minim2.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md1.tpr" -maxwarn 1
    gmx mdrun -deffnm "${prefix}_bilayer_md1"
    gmx grompp -f md2.mdp -c "${prefix}_bilayer_md1.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md2.tpr" -maxwarn 1
    gmx mdrun -deffnm "${prefix}_bilayer_md2 -nb gpu -bonded gpu -pme gpu"
    gmx grompp -f md3.mdp -c "${prefix}_bilayer_md2.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md3.tpr" -maxwarn 1
    gmx mdrun -deffnm "${prefix}_bilayer_md3 -nb gpu -bonded gpu -pme gpu"
    gmx grompp -f md4.mdp -c "${prefix}_bilayer_md3.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md4.tpr" -maxwarn 1
done

