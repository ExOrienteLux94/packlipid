import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from datetime import datetime
import os
import shutil

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



def find_extreme_coords(pdb_file):
    """
    Find the minimum and maximum x, y, z coordinates in a PDB file and the atom indices where they occur.

    Returns:
        extremes (np.ndarray): 3x2 array of min and max for x, y, z
        indices (dict): indices of atoms with min/max x, y, z
    """
    coords = []
    atom_indices = []

    with open(pdb_file, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ndx = int(line[6:11])
                    coords.append([x, y, z])
                    atom_indices.append(ndx)
                except ValueError:
                    continue

    coords = np.array(coords)
    if coords.size == 0:
        raise ValueError("No atom coordinates found in the PDB file.")

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

def main():


    parser = argparse.ArgumentParser(prog='Packlipid.py', description='Python program for writing packmol input scripts for complex lipid bilayer compostions with possibility to include imbedded proteins. Written by Maisa Vuorte\n')
    arg_group_general = parser.add_argument_group('General')
    arg_group_utility = parser.add_argument_group('Utility')
    arg_group_protein = parser.add_argument_group('Protein')
    
    arg_group_general.add_argument("-out", "--out", type=str, default='packmol_input.inp', help='name of packmol output file (default: %(default)s)')
    arg_group_general.add_argument("-comp", "--comp", type=str, required=True, help='name of lipid composition .csv file (default: %(default)s)')
    arg_group_protein.add_argument("-prot", "--prot",type=str, default = 'AF-P11166-F1-model_v4.pdb', help='name of protein struture file (default: %(default)s)')
    arg_group_general.add_argument("-sol", "--sol", type=str, default='./pdb/SOL.pdb', help='name of solvent molecule strcture file (default: %(default)s)')
    arg_group_general.add_argument("-struct", "--struct", choices=['pdb', 'gro'], default="pdb", help='type of structure file to read gro or pdb (default: %(default)s)')
    arg_group_general.add_argument("-pref", "--pref", type=str, default='out', help='prefix of outut strcture file (default: %(default)s)')
    arg_group_general.add_argument("-delim", "--delim", type =str, default=';', help='delimiter used in .csv file (default: %(default)s)')


    # solvation and utility arguments
    arg_group_utility.add_argument("-solvate", "--solvate", type=bool, default=False, help='solvate created strcutre with solvent True/False (default: %(default)s)')
    arg_group_utility.add_argument("-dsol", "--dsol", type=float, default=20, help='height of water layer above and below the bilayer [Å] (default: %(default)s)')
    arg_group_utility.add_argument("-neutral", "--neutral", type=bool, default=False, help='add ions (Na+/Cl-) to neutralize system net charge to 0 True/False (default: %(default)s)')

    # protein insert arguments
    arg_group_protein.add_argument("-proteins", "--proteins", type=bool, default=False, help='inbed protein structures in bilayer True/False (default: %(default)s)')
    arg_group_protein.add_argument("-nprotein", "--nprotein", type=int, help='number of protein structures to inbed in bilayer (default: %(default)s)')
    arg_group_protein.add_argument("-method", "--method", choices=['grid', 'polygon', 'random'], default='grid', help='method for distributing protien structures in bilayer, choises: %(choices)s (default: %(default)s)')

    args = parser.parse_args()

    cwd = os.getcwd()

    delim = args.delim

    struct_out_type = args.struct # type of output file, should add check to see if particle number exceeds 9999 (limit of pdb format)
    struct_out_pref = args.pref # output file prefix
    struct_out_file = struct_out_pref + '.' + struct_out_type

    packmol_out = args.out
    comp_in = args.comp
    protein_pdb = args.prot
    #protein_pdb = 'AF-P11166-F1-model_v4.pdb'
    water_pdb = args.sol

    pbc = True
    protein = args.proteins # embed protein to bilayer
    n_protein = args.nprotein # number of protein copies to include
    protein_distr_mode = args.method # distribution mode 'random', 'grid', or 'polygon'
    plot_protein_dstr = True # plot protein COMs

    solvate = args.solvate # solvate system with water beads
    neutral = args.neutral # include neutralizing ions

    write_topol = True

    run_packmol = False

    d_water = args.dsol # water layer above and below membrane [Å]


    try:
        # Create output directory
        timestamp = datetime.now().strftime("%M-%H-%d-%Y")
        timestamp = datetime.now().strftime("%M-%H-%d-%Y")
        base_dirname = f"run_{timestamp}"
        dirname = base_dirname

        i = 1
        while os.path.exists(dirname):
            dirname = f"{base_dirname}_{i}"
            i += 1

        os.makedirs(dirname)
        print(f"Created output directory: {os.path.join(cwd, dirname)}")

        # copy input files to output directory
        files_to_copy = [comp_in, protein_pdb]

        for src_file in files_to_copy:
            if os.path.exists(src_file):
                dst_file = os.path.join(dirname, os.path.basename(src_file))
                shutil.copy(src_file, dst_file)
                print(f"Copied {src_file} to {dst_file}")
            else:
                print(f"File not found: {src_file}")
        os.chdir(os.path.join(cwd, dirname))
        cwd = os.getcwd()

        # Load composition CSV
        composition = pd.read_csv(comp_in, delimiter=delim, header=None)
        composition.columns = ['lipid', 'upper', 'lower']
        # Load areas CSV
        areas = pd.read_csv('areas.csv', delimiter=delim)

        # Merge on lipid name
        merged = composition.merge(areas, left_on='lipid', right_on='lipid')

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

    except IOError:
        print("One or more of the provided files does not exists or unable to read file contents!")
        quit

    try:
        if protein:
            prot_coords, prot_ndxs = find_extreme_coords(protein_pdb)
            diff_xy = np.abs(prot_coords[0:-1,1] - prot_coords[0:-1,0])
            r_cylider = np.max(diff_xy) * 0.525
            h_cylider = np.abs(prot_coords[-1,1] - prot_coords[1,0]) * 1.05
            prot_area = diff_xy[0] * diff_xy[1] * np.pi * n_protein
            xy_length = (avg_leaflet_area + prot_area) ** 0.5
            print(f"After adding {n_protein} proteins, Box XY dimensions are: {xy_length:.2f} x {xy_length:.2f} Å")

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
        
        # Fixed box Z (change to user query)
        z = d_water * 2 + 21
        z_center = z / 2
        

        with open(packmol_out, 'w') as f:
            f.write("tolerance 3.0\n")
            f.write(f"filetype {struct_out_type}\n")

            if pbc:
                f.write(f"pbc {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
                f.write("add_box_sides 1.0\n")

            f.write(f"output {struct_out_file}\n\n")

            if protein:
                for i in range(n_protein):
                    z_cylinder_start = z_center - h_cylider / 2
                    z_cylinder_end = z_center + h_cylider / 2
                    z_max = prot_ndxs["z_max"]
                    z_min = prot_ndxs["z_min"]
                    f.write(f"structure {protein_pdb}\n")
                    f.write(f"  number 1\n")
                    f.write(f"  inside cylinder {prot_coms[i,0]:.2f} {prot_coms[i,1]:.2f} {z_cylinder_start:.2f} 0. 0. 1. {r_cylider:.2f} {h_cylider:.2f}\n")
                    f.write(f"  atoms {z_max:d}\n")
                    f.write(f"    over plane 0. 0. 1. {z_center+0.5:.2f}\n")
                    f.write("  end atoms\n")
                    f.write(f"  atoms {z_min:d}\n")
                    f.write(f"    below plane 0. 0. 1. {z_center-0.5:.2f}\n")
                    f.write("  end atoms\n")
                    f.write("end structure\n\n")


            for idx, row in merged.iterrows():
                lipid = row['lipid']
                upper = int(row['upper'])
                lower = int(row['lower'])
                pdb_file = f"./pdb/{lipid}-em.pdb"

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
    except SystemExit:
        print("Program execution terminated!")

main()