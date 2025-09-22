#!/bin/bash

#SBATCH --partition=small-g
#SBATCH --account=<project>
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7

module use /appl/local/csc/modulefiles
module load gromacs/2025.2-gpu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

for file in ./*.inp; do
    fname=$(basename "$file" .inp)
    prefix=${fname%%_*}

    if [[ -f "./${prefix}_bilayer_md3.gro" ]]; then
        continue
    fi

    echo "Processing $prefix ..."

    if [[ -f "./${prefix}_bilayer_md2.gro" ]]; then
    	echo -e "del 1 - 20\n r SOL \n r ${prefix} \n name 1 SOL \n name 2 MEMB \n q\n" | srun gmx_mpi make_ndx -f "${prefix}_bilayer.gro" -o "${prefix}_bilayer.ndx"
    	srun gmx_mpi grompp -f md3.mdp -c "${prefix}_bilayer_md2.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md3.tpr" -maxwarn 1
    	srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_md3" -nb gpu -bonded gpu -pme gpu
    	srun gmx_mpi grompp -f md4.mdp -c "${prefix}_bilayer_md3.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md4.tpr" -maxwarn 1
        continue
    fi
	
    if [[ -f "./${prefix}_bilayer_md1.gro" ]]; then
    	echo -e "del 1 - 20\n r SOL \n r ${prefix} \n 0 & ! 1 & ! 2 \n name 1 SOL \n name 2 MEMB \n name 3 Protein \n q\n" | srun gmx_mpi make_ndx -f "${prefix}_bilayer.gro" -o "${prefix}_bilayer.ndx"
	srun gmx_mpi grompp -f md2.mdp -c "${prefix}_bilayer_md1.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md2.tpr" -maxwarn 1
    	srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_md2" -nb gpu -bonded gpu -pme gpu
    	srun gmx_mpi grompp -f md3.mdp -c "${prefix}_bilayer_md2.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md3.tpr" -maxwarn 1
    	srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_md3" -nb gpu -bonded gpu -pme gpu
	srun gmx_mpi grompp -f md4.mdp -c "${prefix}_bilayer_md3.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md4.tpr" -maxwarn 1
        continue
    fi

    srun gmx_mpi editconf -f "${prefix}.pdb" -o "${prefix}_bilayer.gro"
    echo -e "del 1 - 20\n r SOL \n r ${prefix} \n 0 & ! 1 & ! 2 \n name 1 SOL \n name 2 MEMB \n name 3 Protein \n q\n" | srun gmx_mpi make_ndx -f "${prefix}_bilayer.gro" -o "${prefix}_bilayer.ndx"
    srun gmx_mpi grompp -f minim0.mdp -c "${prefix}_bilayer.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_minim0.tpr" -maxwarn 1
    srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_minim0" -nb gpu -bonded gpu -pme gpu
    srun gmx_mpi grompp -f minim1.mdp -c "${prefix}_bilayer_minim0.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_minim1.tpr" -maxwarn 1
    srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_minim1" -nb gpu -bonded gpu -pme gpu
    srun gmx_mpi grompp -f minim2.mdp -c "${prefix}_bilayer_minim1.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_minim2.tpr" -maxwarn 1
    srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_minim2" -nb gpu -bonded gpu -pme gpu
    srun gmx_mpi grompp -f md1.mdp -c "${prefix}_bilayer_minim2.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md1.tpr" -maxwarn 1
    srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_md1" -nb gpu -bonded gpu -pme gpu
    srun gmx_mpi grompp -f md2.mdp -c "${prefix}_bilayer_md1.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md2.tpr" -maxwarn 1
    srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_md2" -nb gpu -bonded gpu -pme gpu
    srun gmx_mpi grompp -f md3.mdp -c "${prefix}_bilayer_md2.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md3.tpr" -maxwarn 1
    srun gmx_mpi mdrun -deffnm "${prefix}_bilayer_md3" -nb gpu -bonded gpu -pme gpu
    srun gmx_mpi grompp -f md4.mdp -c "${prefix}_bilayer_md3.gro" -p "${prefix}_bilayer.top" -n "${prefix}_bilayer.ndx" -r "${prefix}_bilayer.gro" -o "${prefix}_bilayer_md4.tpr" -maxwarn 1
done
