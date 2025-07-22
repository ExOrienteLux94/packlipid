import pandas as pd

# Load composition CSV
composition = pd.read_csv('composition.csv', delimiter=';', header=None)
composition.columns = ['lipid', 'upper', 'lower']

# Load areas CSV
areas = pd.read_csv('areas.csv', delimiter=';')

# Merge on lipid name
merged = composition.merge(areas, left_on='lipid', right_on='Lipid')

# Calculate total leaflet areas
merged['upper_area'] = merged['upper'] * merged['Area']
merged['lower_area'] = merged['lower'] * merged['Area']

# Sum areas
total_upper_area = merged['upper_area'].sum()
total_lower_area = merged['lower_area'].sum()

# Average box xy dimensions (square box, can be changed to allow for general rectangular boxes)
avg_leaflet_area = (total_upper_area + total_lower_area) / 2
xy_length = avg_leaflet_area ** 0.5

# Calculate net charge for adding counter ions
total_charge = (merged['upper'] + merged['lower']) * merged['Charge']
total_charge = total_charge.sum()

# print out 
print(f"Total upper leaflet area: {total_upper_area:.2f} Å²")
print(f"Total lower leaflet area: {total_lower_area:.2f} Å²")
print(f"Box XY dimensions: {xy_length:.2f} x {xy_length:.2f} Å")
print(f"Total bilayer charge: {total_charge:.2f} e")

# Fixed box Z (change to user query)
z = 100
z_center = z / 2

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

with open('packmol_input.inp', 'w') as f:
    f.write("tolerance 2.0\n")
    f.write("filetype pdb\n")
    f.write(f"pbc {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
    f.write("add_box_sides 1.0\n")
    f.write("output bilayer.pdb\n\n")

    for idx, row in merged.iterrows():
        lipid = row['lipid']
        upper = int(row['upper'])
        lower = int(row['lower'])
        pdb_file = f"{lipid}-em.pdb"

        info = find_tail_atoms(pdb_file)
        first = info['first_atom']
        last_A = info['last_A']
        last_B = info['last_B']
        last = info['last_atom']

        # Upper leaflet
        if upper > 0:
            f.write(f"structure {pdb_file}\n")
            f.write(f"  number {upper}\n")
            f.write(f"  inside box 0. 0. {z_center:.2f} {xy_length:.2f} {xy_length:.2f} {z:.2f}\n")
            f.write(f"  restrain {first} {first} z 60.0 65.0\n") # change limits to not hard coded

            if last_A:
                f.write(f"  restrain {last_A} {last_A} z 0.0 55.0\n")
            if last_B:
                f.write(f"  restrain {last_B} {last_B} z 0.0 55.0\n")
            if not last_A and not last_B:
                f.write(f"  restrain {last} {last} z 0.0 55.0\n")

            f.write("end structure\n\n")

        # Lower leaflet
        if lower > 0:
            f.write(f"structure {pdb_file}\n")
            f.write(f"  number {lower}\n")
            f.write(f"  inside box 0. 0. 0. {xy_length:.2f} {xy_length:.2f} {z_center:.2f}\n")
            f.write(f"  restrain {first} {first} z 35.0 40.0\n")

            if last_A:
                f.write(f"  restrain {last_A} {last_A} z 45.0 100.0\n")
            if last_B:
                f.write(f"  restrain {last_B} {last_B} z 45.0 100.0\n")
            if not last_A and not last_B:
                f.write(f"  restrain {last} {last} z 45.0 100.0\n")

            f.write("end structure\n\n")

print("PACKMOL input generated with updated box size and bilayer charge calculated.")