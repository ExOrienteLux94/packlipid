import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from datetime import datetime
import os
import glob
import shutil


def find_tail_atoms(pdb_path):
    tail_A = []
    tail_B = []
    P_atoms = []
    last_atom = None

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()
                last_atom = atom_num
                if 'P' in atom_name:
                    P_atoms.append(atom_num)
                elif 'A' in atom_name:
                    tail_A.append(atom_num)
                elif 'B' in atom_name:
                    tail_B.append(atom_num)

    return {
        'first_atom': min(P_atoms) if P_atoms else None,
        'last_A': tail_A[-1] if tail_A else None,
        'last_B': tail_B[-1] if tail_B else None,
        'last_atom': last_atom
    }

def write_topology(top_file, lipid, system_name='single lipid', lipid_itp=None):
    """
    Write a Martini .top file given a dataframe of lipid composition (upper/lower leaflet)
    and net system charge.

    Args:
        output_file (str): Path to write .top file
        system_name (str): Name of the system
        df (pd.DataFrame): DataFrame with columns ['lipid', 'upper', 'lower']
        total_charge (int): Net charge of system before neutralization
    """

    with open(top_file, "w") as f:
        # Header
        f.write("; Generated Martini topology\n")
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0.itp"\n\n')
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0_ffbonded_v2.itp"\n')
        f.write('#include "/worktmp/scd287/daily/packlipid/lipid_itps/'+ lipid_itp +'"\n')
        f.write('\n')
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0_ions_v1.itp"\n\n')
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0_solvents_v1.itp"\n\n')

        # System section
        f.write("[ system ]\n")
        f.write(f"{system_name}\n\n")

        # Molecules section
        f.write("[ molecules ]\n")
        f.write("; Compound        #mols\n")

        f.write(f"{lipid:<15}{1}\n")

    print(f"Topology written to {top_file}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python process_gro.py lipid_name gro_file_prefix em.mdp")
        sys.exit(1)

    pdb_file = sys.argv[2] + '_bb.pdb'
    gro_file = sys.argv[2] + '_bb.gro'
    mdp_file = sys.argv[3]
    top_file = sys.argv[2] + '.top'
    lipid = sys.argv[1]
    itp = sys.argv[1] + '.itp'

    ndx = find_tail_atoms(pdb_file)
    ndx_file = sys.argv[1]+'.ndx'

    write_topology(top_file, lipid, lipid_itp=itp)

    # create ndx file
    if ndx['last_A'] != None and ndx['last_B'] != None:
        cmd = 'echo -e \" del 1-10 \n 0 & a '+str(ndx['first_atom']) +' \n 0 & a '+str(ndx['last_A'])+' \n 0 & a '+str(ndx['last_B'])+' \n name 1 P \n name 2 tA \n name 3 tB \n q\" | gmx make_ndx -f '+gro_file+' -o '+ndx_file
        nr = 3
    else: 
        cmd = 'echo -e \" del 1-10 \n 0 & a '+str(ndx['first_atom']) +' \n 0 & a '+str(ndx['last_atom'])+'name 1 P \n name 2 tA \n q\" | gmx make_ndx -f '+gro_file+' -o '+ndx_file
        nr = 2
    os.system(cmd)

    if nr == 2:
        mdp_file = 'pull1.mdp'

    tpr_file = sys.argv[1]+'_pull.tpr'
    cmd = 'gmx grompp -f '+mdp_file+' -c '+gro_file+ ' -n '+ndx_file+' -p '+top_file+' -o '+tpr_file

    os.system(cmd)
    cmd = 'gmx mdrun -v -s '+tpr_file
    os.system(cmd)

    gro_out= lipid + '_c.gro'
    pdb_out = lipid + '_c.pdb'
    cmd = 'gmx editconf -f confout.gro -c -box 0 0 0 -o '+gro_out
    os.system(cmd)
    cmd = 'gmx editconf -f '+gro_out+' -o '+ pdb_out
    os.system(cmd)
    os.system(r'rm -rf *#')

main()




