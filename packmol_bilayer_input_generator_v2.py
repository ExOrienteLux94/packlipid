import pandas as pd
import numpy as np


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


def find_extreme_coords(pdb_path, n=1):
    """
    Parses a PDB file and returns a numpy array with the averages of the n maximum and minimum X, Y, Z coordinates.

    Shape: (3, 2) -> [[x_min, x_max], [y_min, y_max], [z_min, z_max]]
    """
    x_coords, y_coords, z_coords = [], [], []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 54:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)
                except ValueError:
                    continue  # Skip malformed lines

    if not x_coords or not y_coords or not z_coords:
        raise ValueError("No valid ATOM/HETATM coordinates found in file.")

    x_max_avg = np.mean(sorted(x_coords, reverse=True)[:n])
    x_min_avg = np.mean(sorted(x_coords)[:n])

    y_max_avg = np.mean(sorted(y_coords, reverse=True)[:n])
    y_min_avg = np.mean(sorted(y_coords)[:n])

    z_max_avg = np.mean(sorted(z_coords, reverse=True)[:n])
    z_min_avg = np.mean(sorted(z_coords)[:n])

    coords = np.array([
        [x_min_avg, x_max_avg],
        [y_min_avg, y_max_avg],
        [z_min_avg, z_max_avg]
    ])

    return coords

def distribute_coms(x_max, y_max, n_points, max_iter=10000, delta_s = 0.01):
    """
    Distribute n points in a 2D box with periodic boundary conditions
    to maximize the minimum point-to-point distance.

    Args:
        x_range (float): maximum x coordinate
        y_range (float): maximum y coordinate
        n_points (int): number of points
        max_iter (int): number of optimization steps
        delta_s (float): iteration step size (dx and dy)

    Returns:
        numpy.ndarray: Array of shape (n_points, 2) with point coordinates
    """
    # Initialize points randomly
    points = np.zeros((n_points, 2))
    points[:, 0] = np.random.uniform(0, x_max, n_points)
    points[:, 1] = np.random.uniform(0, y_max, n_points)

    def pbc_distance(p1, p2):
        dx = np.abs(p1[0] - p2[0])
        dy = np.abs(p1[1] - p2[1])
        dx = min(dx, box_x - dx)
        dy = min(dy, box_y - dy)
        return np.sqrt(dx ** 2 + dy ** 2)

    # Simple repulsive optimization
    for _ in range(max_iter):
        for i in range(n_points):
            force = np.zeros(2)
            for j in range(n_points):
                if i == j:
                    continue
                dx = points[i][0] - points[j][0]
                dy = points[i][1] - points[j][1]
                dx -= x_max * np.round(dx / x_max)
                dy -= y_max* np.round(dy / y_max)
                r2 = dx ** 2 + dy ** 2 + 1e-12
                force[0] += dx / r2
                force[1] += dy / r2

            # Small step in repulsive direction
            points[i][0] += delta_s * force[0]
            points[i][1] += delta_s * force[1]

            # Apply PBC
            points[i][0] = x_range[0] + (points[i][0] - x_range[0]) % box_x
            points[i][1] = y_range[0] + (points[i][1] - y_range[0]) % box_y

    return points

def main():


    delim = ';'
    packmol_out = 'packmol_input.inp'
    comp_in = 'composition.csv'
    protein_pdb = 'AF-P11166-F1-model_v4.pdb'

    pbc = True
    protein = False
    n_protein = 1

    try:
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

    # Fixed box Z (change to user query)
    z = 100
    z_center = z / 2

    if protein:
        coords = find_extreme_coords(protein_pdb)
        if n_protein == 1:
            points = np.array([[xy_length/2, xy_length/2]])
        else:
            points = distribute_coms(xy_length, n_protein)

    with open(packmol_out, 'w') as f:
        f.write("tolerance 2.0\n")
        f.write("filetype pdb\n")

        if pbc:
            f.write(f"pbc {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
            f.write("add_box_sides 1.0\n")

        f.write("output bilayer.pdb\n\n")

        if protein:
            for i in range(n_protein):
                f.write(f"structure {protein_pdb}\n")
                f.write(f"  number 1\n")


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

    print("PACKMOL input generated with updated box size and bilayer charge calculated.")

main()