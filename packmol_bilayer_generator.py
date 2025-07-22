import pandas as pd

# Load CSV
composition = pd.read_csv('composition.csv', delimiter=';', header=None)
composition.columns = ['lipid', 'upper', 'lower']

x, y, z = 500, 500, 100
z_center = 50

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
    f.write("output bilayer.pdb\n\n")

    for idx, row in composition.iterrows():
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
            f.write(f"  inside box 0. 0. {z_center} {x} {y} {z}\n")

            f.write(f"  restrain {first} {first} z 60.0 65.0\n")

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
            f.write(f"  inside box 0. 0. 0. {x} {y} {z_center}\n")

            f.write(f"  restrain {first} {first} z 35.0 40.0\n")

            if last_A:
                f.write(f"  restrain {last_A} {last_A} z 45.0 100.0\n")
            if last_B:
                f.write(f"  restrain {last_B} {last_B} z 45.0 100.0\n")
            if not last_A and not last_B:
                f.write(f"  restrain {last} {last} z 45.0 100.0\n")

            f.write("end structure\n\n")

print("PACKMOL input script generated as 'packmol_input.inp'")