import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from datetime import datetime
import os
import glob
import shutil



atomic_masses = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
    "P": 30.974, "S": 32.06, "MG": 24.305, "ZN": 65.38,
    "FE": 55.845, "CA": 40.078, "K": 39.098, "NA": 22.990,
    "CL": 35.45
}

def find_tail_atoms(pdb_path):
    tail_A = []
    tail_B = []
    last_atom = None

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()
                last_atom = atom_num

                if 'A' in atom_name:
                    tail_A.append(atom_num)
                elif 'B' in atom_name:
                    tail_B.append(atom_num)

    return {
        'first_atom': 1,
        'last_A': tail_A[-1] if tail_A else None,
        'last_B': tail_B[-1] if tail_B else None,
        'last_atom': last_atom
    }

def average_lipid_length(pdb_file):
    """
    Compute the average distance between the phosphorus atom (P*) 
    and the last tail atoms (*A, *B). 

    If no P* is found, use the first atom in the file.
    If no tail atoms are found, return None.
    """
    atoms = []

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_name, np.array([x, y, z])))

    if not atoms:
        raise ValueError("No atoms found in PDB file.")

    # Find reference atom (P*, else first atom)
    ref_atom = None
    for atom_name, coord in atoms:
        if atom_name.startswith("P"):
            ref_atom = coord
            break
    if ref_atom is None:
        ref_atom = atoms[0][1]  # fallback to first atom

    # Find last tail atoms (*A, *B)
    tail_coords = {}
    for atom_name, coord in reversed(atoms):  # reverse to get last occurrence
        if atom_name.endswith("A") and "A" not in tail_coords:
            tail_coords["A"] = coord
        elif atom_name.endswith("B") and "B" not in tail_coords:
            tail_coords["B"] = coord
        if len(tail_coords) == 2:
            break

    # Compute average distance
    if tail_coords:
        distances = [np.linalg.norm(ref_atom - coord) for coord in tail_coords.values()]
        return float(np.mean(distances))
    else:
        return None


def find_extreme_coords(pdb_file, z_range=None):
    """
    Find the minimum and maximum x, y, z coordinates in a PDB file 
    within an optional z-coordinate range, and return atom indices.

    Args:
        pdb_file (str): Path to PDB file
        z_range (tuple or None): (z_min, z_max). If None, use all atoms.

    Returns:
        extremes (np.ndarray): 3x2 array of min and max for x, y, z
        indices (dict): atom indices of min/max x, y, z
    """
    coords = []
    atom_indices = []

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ndx = int(line[6:11])

                    # Apply z-range filter if provided
                    if z_range is None or (z_range[0] <= z <= z_range[1]):
                        coords.append([x, y, z])
                        atom_indices.append(ndx)
                except ValueError:
                    continue

    coords = np.array(coords)
    if coords.size == 0:
        raise ValueError("No atom coordinates found in the PDB file (or all filtered by z_range).")

    extremes = np.zeros((3, 2))
    indices = {'x_min': None, 'x_max': None, 'y_min': None, 'y_max': None, 'z_min': None, 'z_max': None}

    for i, axis in enumerate(['x', 'y', 'z']):
        min_val = np.min(coords[:, i])
        max_val = np.max(coords[:, i])
        extremes[i, 0] = min_val
        extremes[i, 1] = max_val

        indices[f'{axis}_min'] = atom_indices[np.argmin(coords[:, i])]
        indices[f'{axis}_max'] = atom_indices[np.argmax(coords[:, i])]

    return extremes, indices

def distribute_coms(x_range, y_range, n_points, max_iter=10000, delta_s = 1):
    """
    Distribute n points in a 2D box with periodic boundary conditions
    to maximize the minimum point-to-point distance using gradient-based optimization.

    Args:
        x_range (np.ndarray): (xmin, xmax)
        y_range (np.ndarray): (ymin, ymax)
        n_points (int): number of points
        max_iter (int): number of optimization steps
        lr (float): learning rate for optimization

    Returns:
        numpy.ndarray: Array of shape (n_points, 2) with point coordinates
    """
    # Initialize points randomly
    # Initialize points randomly
    points = np.zeros((n_points, 2))
    points[:, 0] = np.random.uniform(x_range[0], x_range[1], n_points)
    points[:, 1] = np.random.uniform(y_range[0], y_range[1], n_points)

    box_x = x_range[1] - x_range[0]
    box_y = y_range[1] - y_range[0]

    def pbc_vector(p1, p2):
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]
        dx -= box_x * np.round(dx / box_x)
        dy -= box_y * np.round(dy / box_y)
        return np.array([dx, dy])

    # Gradient-based repulsion optimization
    for _ in range(max_iter):
        gradients = np.zeros_like(points)
        for i in range(n_points):
            for j in range(n_points):
                if i == j:
                    continue
                delta = pbc_vector(points[i], points[j])
                dist_sq = np.sum(delta ** 2) + 1e-12
                gradients[i] += delta / dist_sq

        # Gradient descent step
        points += delta_s * gradients

        # Apply PBC
        points[:, 0] = x_range[0] + (points[:, 0] - x_range[0]) % box_x
        points[:, 1] = y_range[0] + (points[:, 1] - y_range[0]) % box_y

    return points

def grid_points_pbc(x_range, y_range, n_points):
    """
    Place n points in a 2D periodic box with maximum equal spacing using a regular grid.

    Args:
        x_range (tuple): (xmin, xmax)
        y_range (tuple): (ymin, ymax)
        n_points (int): Number of points

    Returns:
        np.ndarray: (n_points, 2) array of point coordinates
    """
    box_x = x_range[1] - x_range[0]
    box_y = y_range[1] - y_range[0]

    # Estimate grid size (rows x cols) as close to square as possible
    n_cols = int(np.ceil(np.sqrt(n_points)))
    n_rows = int(np.ceil(n_points / n_cols))

    # Adjust spacing based on grid dimensions
    spacing_x = box_x / n_cols
    spacing_y = box_y / n_rows

    # Generate points
    points = []
    for i in range(n_rows):
        for j in range(n_cols):
            if len(points) >= n_points:
                break
            x = x_range[0] + (j + 0.5) * spacing_x
            y = y_range[0] + (i + 0.5) * spacing_y
            points.append((x, y))

    return np.array(points)

def regular_polygon_pbc(x_range, y_range, n_points):
    """
    Place n points in a regular polygon arrangement inside a periodic box.
    Maximizes point-to-point distances while avoiding PBC overlap.

    Args:
        x_range (tuple): (xmin, xmax)
        y_range (tuple): (ymin, ymax)
        n_points (int): number of points (>=3 recommended)

    Returns:
        np.ndarray: (n_points, 2) array of coordinates
    """
    box_x = x_range[1] - x_range[0]
    box_y = y_range[1] - y_range[0]

    # Center of box
    center_x = (x_range[0] + x_range[1]) / 2
    center_y = (y_range[0] + y_range[1]) / 2

    # Maximum safe radius without overlapping across PBCs
    max_radius = min(box_x, box_y) / 3

    # Angle between points
    angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)

    # Coordinates of points
    x_coords = center_x + max_radius * np.cos(angles)
    y_coords = center_y + max_radius * np.sin(angles)

    points = np.stack([x_coords, y_coords], axis=1)

    return points



def evaluate_distribution(points, x_range, y_range):
    """
    Evaluates and plots the distribution of points.
    """
    n_points = len(points)
    box_x = x_range[1] - x_range[0]
    box_y = y_range[1] - y_range[0]

    dists = []
    for i in range(n_points):
        for j in range(i + 1, n_points):
            dx = points[i][0] - points[j][0]
            dy = points[i][1] - points[j][1]
            dx -= box_x * np.round(dx / box_x)
            dy -= box_y * np.round(dy / box_y)
            dist = np.sqrt(dx ** 2 + dy ** 2)
            dists.append(dist)

    min_dist = np.min(dists)
    avg_dist = np.mean(dists)

    print(f"Minimum distance between protein COMs: {min_dist:.3f}")
    print(f"Average distance between protein COMs: {avg_dist:.3f}")

    plt.figure(figsize=(6, 6))
    plt.scatter(points[:, 0], points[:, 1], s=30, color='blue')
    plt.xlim(x_range)
    plt.ylim(y_range)
    plt.gca().set_aspect('equal')
    plt.title("Point Distribution with PBC (Gradient Optimized)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.show()

def get_mass(element):
    """Return atomic mass for element, default 12.0 if unknown"""
    return atomic_masses.get(element.upper(), 12.0)

def center_pdb(input_pdb, output_pdb, box_size):
    """
    Center a protein PDB structure by COM and set box size.

    Args:
        input_pdb (str): input PDB filename
        output_pdb (str): output PDB filename
        box_size (tuple): (x, y, z) box dimensions
    """
    coords = []
    masses = []
    lines = []


    # Parse coordinates + elements
    with open(input_pdb, "r") as f:
        for line in f:
            lines.append(line)  # store ALL lines
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

                # Try to infer element
                element = line[76:78].strip()
                if not element:
                    element = line[12:14].strip()[0]
                masses.append(get_mass(element))

    coords = np.array(coords)
    masses = np.array(masses)

    if coords.size == 0:
        raise ValueError("No atoms found in PDB file.")

    # Compute center of mass
    com = np.average(coords, axis=0, weights=masses)

    # Target center of box
    box_center = np.array([box_size[0]/2, box_size[1]/2, box_size[2]/2])

    # Shift vector
    shift = box_center - com

    # Write new PDB with shifted coordinates
    with open(output_pdb, "w") as f:
        # Always insert new CRYST1 as first line
        cryst1 = f"CRYST1{box_size[0]:9.3f}{box_size[1]:9.3f}{box_size[2]:9.3f}" \
                 f"{90.00:7.2f}{90.00:7.2f}{90.00:7.2f} P 1           1\n"
        f.write(cryst1)

        atom_idx = 0
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                x = coords[atom_idx, 0] + shift[0]
                y = coords[atom_idx, 1] + shift[1]
                z = coords[atom_idx, 2] + shift[2]
                new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                f.write(new_line)
                atom_idx += 1
            elif not line.startswith("CRYST1"):  # skip old CRYST1 if present
                f.write(line)

    #print(f"Protein centered in box {box_size} (COM method), written to {output_pdb}")

def get_molnames(itp_file, target=None):
    molnames = []
    with open(itp_file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines, start=1):
        if line.strip().startswith("; molname"):
            # The next non-empty, non-comment line is the molecule name
            j = i
            while j < len(lines):
                candidate = lines[j].strip()
                if candidate and not candidate.startswith(";") and not candidate.startswith("["):
                    name = candidate.split()[0]
                    if target is None or name == target:
                        molnames.append((j+1, name))  # store line number and name
                    break
                j += 1
    return molnames

def write_topology(top_file, df, total_charge, system_name='bilayer system', lipid_itp=None, protein_itp = None, n_prot = 0, dimer_bool = False, protein_bool = False):
    """
    Write a Martini .top file given a dataframe of lipid composition (upper/lower leaflet)
    and net system charge.

    Args:
        output_file (str): Path to write .top file
        system_name (str): Name of the system
        df (pd.DataFrame): DataFrame with columns ['lipid', 'upper', 'lower']
        total_charge (int): Net charge of system before neutralization
    """
    # Check required columns
    if not {"lipid", "upper", "lower"}.issubset(df.columns):
        raise ValueError("DataFrame must contain 'lipid', 'upper', and 'lower' columns")

    # Determine ions required for neutralization
    ions = {}
    if total_charge > 0:
        ions["CL-"] = total_charge
    elif total_charge < 0:
        ions["NA+"] = -total_charge

    with open(top_file, "w") as f:
        # Header
        f.write("; Generated Martini topology\n")
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0.itp"\n\n')
        for litp in lipid_itp:
            f.write('#include "/worktmp/scd287/daily/martini3/lipids/'+ litp +'"\n')
        f.write('\n')
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0_ions_v1.itp"\n\n')
        f.write('#include "/worktmp/scd287/daily/martini3/lipids/martini_v3.0.0_solvents_v1.itp"\n\n')

        if protein_bool:
            f.write('#include "/worktmp/scd287/daily/martini3/martini_v3.0.0_proteins/force_fields/martini3001/general.ff"\n')
            f.write('#include "/worktmp/scd287/daily/martini3/martini_v3.0.0_proteins/force_fields/martini3001/aminoacids.ff"\n')
            f.write('#include "/worktmp/scd287/daily/martini3/martini_v3.0.0_proteins/force_fields/martini3001/modifications.ff"\n\n')
            for prot in protein_itp:
                f.write('#include "' + prot +'"\n')
            f.write('\n')
        else:
            f.write('\n')

        # System section
        f.write("[ system ]\n")
        f.write(f"{system_name}\n\n")

        # Molecules section
        f.write("[ molecules ]\n")
        f.write("; Compound        #mols\n")
        if protein_bool == True:
            for itp in protein_itp:
                mol_name = get_molname_from_itp(itp)    
                if dimer_bool:
                    f.write(f"{mol_name:<15}{n_prot}\n")
                else:
                    f.write(f"{mol_name:<15}{n_prot}\n")
        for _, row in df.iterrows():
            total_count = int(row["upper"]) + int(row["lower"])
            f.write(f"{row['lipid']:<15}{total_count}\n")
        for ion, count in ions.items():
            f.write(f"{ion:<15}{count}\n")

    print(f"Topology written to {top_file}")
    print(f"Net system charge before neutralization: {total_charge}")
    if ions:
        print(f"Added ions: {ions}")
    else:
        print("System already neutral.")


def get_molname_from_itp(itp_file):
    """
    DEPRECATED!!!!
    Extract the molecule name from a .itp file (from the [ moleculetype ] section).
    """
    with open(itp_file, "r") as f:
        lines = f.readlines()

    mol_name = None
    in_section = False
    for line in lines:
        line = line.strip()
        if line.startswith("[ moleculetype ]"):
            in_section = True
            continue
        if in_section:
            if line.startswith(";") or not line:  # skip comments/empty
                continue
            mol_name = line.split()[0]  # molecule name is first token
            break

    return mol_name

def net_charge_from_itp(itp_file):
    """
    Parse a GROMACS .itp file and compute net molecular charge
    from the [ atoms ] section.
    """
    in_atoms = False
    total_charge = 0.0

    with open(itp_file, "r") as f:
        for line in f:
            line = line.strip()
            # Section handling
            if line.startswith("["):
                if line.lower().startswith("[ atoms ]") or line.lower().startswith("[atoms]"):
                    in_atoms = True
                else:
                    in_atoms = False
                continue

            # Process lines inside [ atoms ]
            if in_atoms and line and not line.startswith(";"):
                parts = line.split()
                try:
                    charge = float(parts[6])  # 7th column is charge
                    total_charge += charge
                except (IndexError, ValueError):
                    continue

    return total_charge

def main():


    parser = argparse.ArgumentParser(prog='Packlipid.py', description='Python program for writing packmol input scripts for complex lipid bilayer compostions with possibility to include imbedded proteins. Written by Maisa Vuorte\n')
    arg_group_general = parser.add_argument_group('General')
    arg_group_utility = parser.add_argument_group('Utility')
    arg_group_protein = parser.add_argument_group('Protein')
    
    arg_group_general.add_argument("-out", "--out", type=str, default='packmol_input.inp', help='name of packmol output file (default: %(default)s)')
    arg_group_general.add_argument("-top", "--top", type=str, default='topol.top', help='name of gromacs topology .top output file (default: %(default)s)')
    arg_group_general.add_argument("-comp", "--comp", type=str, default=None, help='name of lipid composition .csv file (default: %(default)s)')
    arg_group_general.add_argument("-pitp", "--pitp", type=str, nargs='*', default=None, help='name(s) of protein gromacs topology .itp files (default: %(default)s)')
    arg_group_general.add_argument("-sol", "--sol", type=str, default=None, help='name of solvent molecule strcture file (default: %(default)s)')
    arg_group_general.add_argument("-struct", "--struct", choices=['pdb', 'gro'], default="pdb", help='type of structure file to read gro or pdb (default: %(default)s)')
    arg_group_general.add_argument("-pref", "--pref", type=str, default='out', help='prefix of output structure file (default: %(default)s)')
    arg_group_general.add_argument("-delim", "--delim", type =str, default=';', help='delimiter used in .csv file (default: %(default)s)')
    arg_group_general.add_argument("-delim", "--delim", type =str, default=None, help='delimiter used in .csv file (default: %(default)s)')
    arg_group_general.add_argument("-ff", "--ff", type=str, choices=['martini2', 'martini3'], default='martini3', help=' (default: %(default)s)')


    # solvation and utility arguments
    arg_group_utility.add_argument("-lipids", "--lipids", nargs='*', type = str, default=None, help='Name(s) of lipids constituting the membrane when not reading composition from .csv, order must match --nupper and --nlower (default: %(default)s)')
    arg_group_utility.add_argument("-nupper", "--nupper", nargs='*', type = int, default=None, help='Number of lipids constituting the membrane when not reading composition from .csv, order must match --lipids (default: %(default)s)')
    arg_group_utility.add_argument("-nlower", "--nlower", nargs='*', type = int, default=None, help='Number of lipids constituting the membrane when not reading composition from .csv, order must match --lipids (default: %(default)s)')
    arg_group_utility.add_argument("-solvate", "--solvate", type=str, default='no', help='solvate created strcutre with solvent yes/no (default: %(default)s)')
    arg_group_utility.add_argument("-dsol", "--dsol", type=float, default=20, help='height of water layer above and below the bilayer [Å] (default: %(default)s)')
    arg_group_utility.add_argument("-neutral", "--neutral", type=str, default='no', help='add ions (Na+/Cl-) to neutralize system net charge to 0 yes/no (default: %(default)s)')

    # protein insert arguments
    arg_group_protein.add_argument("-prot", "--prot",type=str, default = None, help='name of protein struture file (default: %(default)s)')
    arg_group_protein.add_argument("-proteins", "--proteins", type=str, default='no', help='inbed protein structures in bilayer yes/no (default: %(default)s)')
    arg_group_protein.add_argument("-dimer", "--dimer", type=str, default='no', help='protein .pdb file contains a dimer yes/no (default: %(default)s)')
    arg_group_protein.add_argument("-orient", "--orient", type=str, default='no', help='reorient protein structure with respect to bilayer yes/no (default: %(default)s)')
    arg_group_protein.add_argument("-nprotein", "--nprotein", type=int, help='number of protein structures to inbed in bilayer (default: %(default)s)')
    arg_group_protein.add_argument("-method", "--method", choices=['replace', 'scale'], default='replace', help='method for inserting protein structure in bilayer, options: %(choices)s (default: %(default)s)')
    arg_group_protein.add_argument("-dstr", "--dstr", choices=['grid', 'polygon', 'random'], default='grid', help='method for distributing protein structures in bilayer, options: %(choices)s (default: %(default)s)')

    args = parser.parse_args()

    test = True
    #
    cwd = os.getcwd()
    packlipids_dir = os.path.dirname(os.path.realpath(sys.argv[0])) # directory where packlipids is installed

    comp_from_csv = True
    if args.lipids != None and args.nupper != None and args.nlower != None:
        comp_from_csv = False

    delim = args.delim
    ffield = args.ff

    struct_out_type = args.struct # type of output file, should add check to see if particle number exceeds 9999 (limit of pdb format)
    struct_out_pref = args.pref # output file prefix
    struct_out_file = struct_out_pref + '.' + struct_out_type

    packmol_out = args.out
    comp_in = args.comp
    protein_pdb = args.prot
    #protein_pdb = 'AF-P11166-F1-model_v4.pdb'
    water_pdb = args.sol

    pbc = True
    if args.proteins == 'yes':
        protein = True  # embed protein to bilayer
    else:
        protein = False

    n_protein = args.nprotein # number of protein copies to include
    protein_distr_mode = args.dstr # distribution mode 'random', 'grid', or 'polygon'
    protein_insert_method = args.method # insertion method 'replace' or 'scale'
    if args.dimer == 'yes':
        dimer = True  # protein .pdb is a dimer
    else:
        dimer = False

    plot_protein_dstr = True # plot protein COMs

    if args.solvate == 'yes':
        solvate = True  # solvate system with water beads
    else:
        solvate = False

    if args.neutral == 'yes':
        neutral = True  # include neutralizing ions
    else:
        neutral = False


    write_topol = True
    top_file = args.top
    prot_itps = args.pitp

    rescale_box = True
    center = True # center protein in box

    d_water = args.dsol # water layer above and below membrane [Å]

    d_bilayer = 40 # initial guess for bilayer thickness

    try:
        # Create output directory
        outdir_name = args.outdir

        if outdir_name != None:
            base_dirname = outdir_name

        else:
            timestamp = datetime.now().strftime("%Y-%m-%d-%H%M")
            base_dirname = f"run_{timestamp}"
        
        dirname = base_dirname

        i = 1
        while os.path.exists(dirname) and overwrite == False:
            dirname = f"{base_dirname}_{i}"
            i += 1

        os.makedirs(dirname)
        output_dir = os.path.join(cwd, dirname)
        print(f"Created output directory: {output_dir}")

        # copy input files to output directory
        if comp_in != None:
            files_to_copy = [comp_in, protein_pdb] + prot_itps
        else:
            files_to_copy = [protein_pdb] + prot_itps
        #print(files_to_copy)
        for src_file in files_to_copy:
            if os.path.exists(src_file):
                dst_file = os.path.join(dirname, os.path.basename(src_file))
                shutil.copy(src_file, dst_file)
                print(f"Copied {src_file} to {dst_file}")
            else:
                print(f"File not found: {src_file}")

        if comp_from_csv == True:
            # Load composition CSV
            composition = pd.read_csv(comp_in, delimiter=delim, header=None)
            composition.columns = ['lipid', 'upper', 'lower']
        else:
            if not (len(args.lipids) == len(args.nupper) == len(args.nlower)):
                raise ValueError("Options --lipids, --nupper, and --nlower must have the same number of arguments")
            composition = pd.DataFrame({
                "lipid": args.lipids,
                "upper": args.nupper,
                "lower": args.nlower
            })

        # Load areas CSV
        areas_csv = packlipids_dir + '/data/areas.csv'
        areas = pd.read_csv(areas_csv, delimiter=delim)

        #os.chdir(os.path.join(cwd, dirname))
        #cwd = os.getcwd()

        # Merge on lipid name
        merged = pd.merge(composition, areas, how="left", on="lipid")
        merged["area"]   = merged["area"].fillna(60) # fill out missing areas

        #check for NaN values in charge
        nan_ndx = merged[merged["charge"].isna()].index

        if len(nan_ndx)>0:
            nan_lipids = merged.loc[merged["charge"].isna(), "lipid"].tolist()
            # query lipid net harge from itp file
            charges_from_itp = []
            for ln in nan_lipids:
                itp_dir = packlipids_dir + '/data/ffs/martini3/lipids/top/'
                ln_itp=itp_dir + ln +'.itp'
                if os.path.isfile(ln_itp):
                    charges_from_itp.append(net_charge_from_itp(ln_itp))
                else:
                    # append with default values
                    charges_from_itp.append(0)
            merged.loc[merged["charge"].isna(), "charge"] = charges_from_itp

        # Calculate total leaflet areas
        merged['upper_area'] = merged['upper'] * merged['area']
        merged['lower_area'] = merged['lower'] * merged['area']

        # Sum areas
        total_upper_area = merged['upper_area'].sum()
        total_lower_area = merged['lower_area'].sum()

        # Average box xy dimensions (square box, can be changed to allow for general rectangular boxes)
        avg_leaflet_area = (total_upper_area + total_lower_area) / 2
        xy_length = avg_leaflet_area ** 0.5

        # Calculate net charge for adding counter ions
        total_charge = (merged['upper'] + merged['lower']) * merged['charge']
        total_charge = total_charge.sum()

        # print out 
        print(f"Total upper leaflet area: {total_upper_area:.2f} Å²")
        print(f"Total lower leaflet area: {total_lower_area:.2f} Å²")
        print(f"Box XY dimensions: {xy_length:.2f} x {xy_length:.2f} Å")
        print(f"Total bilayer charge: {total_charge:.2f} e")


        itp_list = []
        if ffield == 'martini3':
            itp_dir = packlipids_dir + '/data/ffs/martini3/lipids/'
            # figure out which lipid its are needed
            for lipid in merged['lipid']:
                suffix = lipid[-2:]
                if suffix == 'PA':
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_PA_v2.itp')
                    itp_list.append(itp_dir+'martini_v3.0.0_etherphospholipids_PA_v2.itp')
                elif suffix == 'PC':
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_PC_v2.itp')
                    itp_list.append(itp_dir+'martini_v3.0.0_etherphospholipids_PC_v2.itp')
                elif suffix == 'PE': 
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_PE_v2.itp')
                    itp_list.append(itp_dir+'martini_v3.0.0_etherphospholipids_PE_v2.itp')
                elif suffix == 'PG': 
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_PG_v2.itp')
                    itp_list.append(itp_dir+'martini_v3.0.0_etherphospholipids_PG_v2.itp')
                elif suffix == 'PS': 
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_PS_v2.itp')
                    itp_list.append(itp_dir+'martini_v3.0.0_etherphospholipids_PS_v2.itp')
                elif suffix == 'PI': 
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_PI_v2.itp')
                elif suffix == 'CL': 
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_CL_v2.itp')
                elif suffix == 'SM': 
                    itp_list.append(itp_dir+'martini_v3.0.0_phospholipids_SM_v2.itp')
                else:
                    lipid_itps = glob.glob(os.path.join(itp_dir, "*.itp"))
                    for fn in lipid_itps:
                        if lipid in get_molnames(fn, target=str(lipid)):
                            itp_list.append(fn)
                            break
            
        #print(lipid_itps)

    except IOError:
        print("One or more of the provided files does not exists or unable to read file contents!")
        quit

    try:
        # determine bilayer width
        lipid_lengths = []
        tot_lipids = 0

        for ndx, row in merged.iterrows():
            lipid = row['lipid']
            nr = int(row['upper']) + int(row['lower'])
            if ffield == 'martini3':
                ll_pdb_file = f"../pdb/martini3/{lipid}.pdb"
            elif ffield =='martini2':
                ll_pdb_file = f"../pdb/martini2/{lipid}-em.pdb"
            ll = average_lipid_length(ll_pdb_file)
            if ll != None:
                lipid_lengths.append(ll*nr)
                tot_lipids += nr

        if tot_lipids != 0:
            d_bilayer = np.sum(lipid_lengths) / tot_lipids * 2

        
        # Determine box z
        z = d_water*2 + d_bilayer

        if protein == True:
            prot_coords, prot_ndxs = find_extreme_coords(protein_pdb)
            diff_xy = np.abs(prot_coords[0:-1,1] - prot_coords[0:-1,0])
            r_cylider = np.max(diff_xy) * 0.525
            h_cylider = np.abs(prot_coords[-1,1] - prot_coords[-1,0]) * 1.05
            if (h_cylider + 40) > z: 
                z = h_cylider +40
            if np.max(diff_xy) > xy_length and rescale_box == True and n_protein == 1:
                # scale mebrane components
                print("!!!!!!!!!!!!!!!!!!!!!!!")  
                print("Initial bilayer dimensions are too small to fit imbedded protein structure with minimum immage convention.")
                print("Rescaling box size to fit protein in membrane.")
                xy_length_old = xy_length 
                xy_length = np.max(diff_xy) * 1.8

            # center protein in box 
            if center:
                centered_out_pdb = output_dir + '/protein_center.pdb'
                center_pdb(protein_pdb, centered_out_pdb, (xy_length, xy_length, z))
            # estimate protein crosssectional area
            cros_prot_coords, cros_prot_ndxs = find_extreme_coords(centered_out_pdb, np.array([z/2 - d_bilayer/2, z/2 + d_bilayer/2]))
            
            cross_diff_xy = np.abs(cros_prot_coords[0:-1,1] - cros_prot_coords[0:-1,0])
            prot_area = (cross_diff_xy [0]*0.5) * (cross_diff_xy [1]*0.5) * np.pi * n_protein

            if protein_insert_method == 'replace' and n_protein > 1:  
                ## NOT COMPLETE YET!!!!!                 
                xy_length = (avg_leaflet_area + prot_area) ** 0.5
                print(f"After adding {n_protein} protein(s), Box XY dimensions are: {xy_length:.2f} x {xy_length:.2f} Å")

            elif protein_insert_method == 'replace' and n_protein == 1:
                # make sure protein fits in box, if not, resize membrane
                if rescale_box == True:
                    area_lipid_new = xy_length ** 2 - prot_area
                    l_scale = (area_lipid_new ** 0.5) / xy_length_old
                    a_scale = int(np.round(l_scale**2, decimals=0))

                    merged['upper'] = merged['upper'] * a_scale
                    merged['lower'] = merged['lower'] * a_scale
                    # Calculate total leaflet areas
                    merged['upper_area'] = merged['upper'] * merged['area']
                    merged['lower_area'] = merged['lower'] * merged['area']
                    # Sum areas
                    total_upper_area = merged['upper_area'].sum()
                    total_lower_area = merged['lower_area'].sum()  
                    # Average box xy dimensions (square box, can be changed to allow for general rectangular boxes)
                    avg_leaflet_area = (total_upper_area + total_lower_area) / 2

                    # Calculate net charge for adding counter ions
                    total_charge = (merged['upper'] + merged['lower']) * merged['charge']
                    total_charge = total_charge.sum()

                    # print out 
                    print("After scaling the new box dimensions are:")
                    print(f"Estimated protein crosssectional area: {prot_area:.2f} Å²")
                    print(f"Total upper leaflet lipid area: {total_upper_area:.2f} Å²")
                    print(f"Total lower leaflet lipid area: {total_lower_area:.2f} Å²")
                    print(f"Box XY dimensions: {xy_length:.2f} x {xy_length:.2f} Å")
                    print(f"Total bilayer charge: {total_charge:.2f} e")
                    print("!!!!!!!!!!!!!!!!!!!!!!!")



            if n_protein == 1:
                prot_coms = np.array([[xy_length/2, xy_length/2]])
            else:
                if protein_distr_mode == 'random':
                    prot_coms = distribute_coms(np.array([0, xy_length]), np.array([0, xy_length]), n_protein)
                elif protein_distr_mode == 'grid':
                    prot_coms = grid_points_pbc(np.array([0, xy_length]), np.array([0, xy_length]), n_protein)
                elif protein_distr_mode == 'polygon':
                    prot_coms = regular_polygon_pbc(np.array([0, xy_length]), np.array([0, xy_length]), n_protein)
                else:
                    print("No valid protein distribution scheme provided as input!")
                    sys.exit()

                if plot_protein_dstr:
                    evaluate_distribution(prot_coms, np.array([0, xy_length]), np.array([0, xy_length]))
        
        z_center = z / 2        

        packmol_out = output_dir + '/'+ packmol_out
        with open(packmol_out, 'w') as f:
            f.write("tolerance 5.0\n")
            f.write(f"filetype {struct_out_type}\n")

            if pbc:
                f.write(f"pbc {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
                f.write("add_box_sides 1.0\n")

            f.write(f"output {struct_out_file}\n\n")

            if (protein == True) and (protein_insert_method == 'replace'):
                if n_protein == 1:
                    #print(prot_coms, diff_xy)
                    f.write(f"structure {centered_out_pdb}\n")
                    f.write(f"  number 1\n")
                    f.write(f"  fixed 0. 0. 0. 0. 0. 0.\n")
                    f.write("end structure\n\n") 
                if n_protein > 1:
                    for i in range(n_protein):
                        z_cylinder_start = z_center - h_cylider / 2
                        z_cylinder_end = z_center + h_cylider / 2
                        z_max = prot_ndxs["z_max"]
                        z_min = prot_ndxs["z_min"]
                        #print(prot_coms, diff_xy)
                        f.write(f"structure {protein_pdb}\n")
                        f.write(f"  number 1\n")
                        f.write(f"  center\n")
                        f.write(f"  fixed {prot_coms[i,0]:.2f} {prot_coms[i,1]:.2f} {z_center:.2f} 0. 0. 0.\n")
                        f.write("end structure\n\n") 

            elif (protein == True) and (protein_insert_method == 'scale'):
                for i in range(n_protein):
                    z_cylinder_start = z_center - h_cylider / 2
                    z_cylinder_end = z_center + h_cylider / 2
                    z_max = prot_ndxs["z_max"]
                    z_min = prot_ndxs["z_min"]
                    #print(prot_coms, diff_xy)
                    f.write(f"structure {protein_pdb}\n")
                    f.write(f"  number 1\n")
                    f.write(f"  inside box {prot_coms[i,0]-(diff_xy[0]*0.5):.2f} {prot_coms[i,1]-(diff_xy[1]*0.5):.2f} {z_cylinder_start:.2f} {prot_coms[i,0]+(diff_xy[0]*0.5):.2f} {prot_coms[i,1]+(diff_xy[1]*0.5):.2f} {z_cylinder_end:.2f}\n")
                    f.write(f"  atoms {z_max:d}\n")
                    f.write(f"    over plane 0. 0. 1. {z_center+10:.2f}\n")
                    f.write("  end atoms\n")
                    f.write(f"  atoms {z_min:d}\n")
                    f.write(f"    below plane 0. 0. 1. {z_center-10:.2f}\n")
                    f.write("  end atoms\n")
                    f.write("end structure\n\n")       

            # build lipid bilayer
            for idx, row in merged.iterrows():
                lipid = row['lipid']
                upper = int(row['upper'])
                lower = int(row['lower'])
                pdb_file = f"../pdb/{lipid}-em.pdb"

                info = find_tail_atoms(pdb_file)
                first = info['first_atom']
                last_A = info['last_A']
                last_B = info['last_B']
                last = info['last_atom']

                # Upper leaflet
                if upper > 0:
                    ul_lim = z_center + 10
                    f.write(f"structure {pdb_file}\n")
                    f.write(f"  number {upper}\n")
                    f.write(f"  inside box 0. 0. {z_center:.2f} {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
                    f.write(f"  atoms {first}\n")
                    f.write(f"    over plane 0. 0. 1. {ul_lim:.2f}\n")
                    f.write("  end atoms\n")
                    
                    if last_A and last_B:
                        f.write(f"  atoms {last_A} {last_B}\n")
                    elif last_A:
                        f.write(f"  atoms {last_A}\n")
                    elif last_B:
                        f.write(f"  atoms {last_B}\n")
                    elif not last_A and not last_B:
                        f.write(f"  atoms {last}\n")
                    f.write(f"    below plane 0. 0. 1. {ul_lim:.2f}\n")
                    f.write("  end atoms\n")
                    f.write("end structure\n\n")

                # Lower leaflet
                if lower > 0:
                    ll_lim = z_center - 10
                    f.write(f"structure {pdb_file}\n")
                    f.write(f"  number {lower}\n")
                    f.write(f"  inside box 0. 0. 0. {xy_length:.2f} {xy_length:.2f} {z_center:.2f}\n")
                    f.write(f"  atoms {first}\n")
                    f.write(f"    below plane 0. 0. 1. {ll_lim:.2f}\n")
                    f.write("  end atoms\n")
                    
                    if last_A and last_B:
                        f.write(f"  atoms {last_A} {last_B}\n")
                    elif last_A:
                        f.write(f"  atoms {last_A}\n")
                    elif last_B:
                        f.write(f"  atoms {last_B}\n")
                    elif not last_A and not last_B:
                        f.write(f"  atoms {last}\n")
                    f.write(f"    over plane 0. 0. 1. {ll_lim:.2f}\n")
                    f.write("  end atoms\n")
                    f.write("end structure\n\n")

            if solvate:
                # assume water density at 1atm, 290 K (55.5 mol/l) and 4 water per MARTINI bead
                rho_w = 8.36  # beads per nm³
                n_water = round((d_water * xy_length * xy_length) * (rho_w / 1000))
                print(f"Solvated system with {2*n_water} water particles.")
                f.write(f"structure {water_pdb}\n")
                f.write(f"  number {n_water}\n")
                f.write(f"  inside box 0. 0. {ul_lim:.2f} {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
                f.write("end structure\n\n")

                f.write(f"structure {water_pdb}\n")
                f.write(f"  number {n_water}\n")
                f.write(f"  inside box 0. 0. 0. {xy_length:.2f} {xy_length:.2f} {ll_lim:.2f}\n")
                f.write("end structure\n\n")

        print("PACKMOL input generated with updated box size and bilayer charge calculated.")

        if write_topol:
            print(lipid_itps)
            top_file = output_dir + '/' + top_file
            write_topology(top_file, df = merged, total_charge=int(total_charge), protein_bool = protein, n_prot= n_protein, dimer_bool= dimer, protein_itp=prot_itps, lipid_itp = lipid_itps)
    except SystemExit:
        print("Program execution terminated!")

main()