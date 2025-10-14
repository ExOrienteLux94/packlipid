import os
import re

def get_atom_names_itp(filename):
    """Extract atom names from [ atoms ] section in .itp file."""
    atoms = []
    in_atoms = False
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                in_atoms = line.startswith('[ atoms ]')
                in_atoms2 = line.startswith('[atoms]')
                continue
            if (in_atoms or in_atoms2) and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 5:
                    atoms.append(parts[4])  # 5th column is atom name
    return atoms

def get_atom_names_gro(filename):
    """Extract atom names from .gro file."""
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        # skip first two lines (title, atom count) and last line (box)
        for line in lines[2:-1]:
            # atom name is in columns 11–15 (fixed width format)
            atom_name = line[10:15].strip()
            atoms.append(atom_name)
    return atoms

# Compare all .itp files
mismatched = []

for itp_file in [f for f in os.listdir('.') if f.endswith('.itp')]:
    base = os.path.splitext(itp_file)[0]
    gro_file = base + '.gro'
    if not os.path.exists(gro_file):
        print(f"⚠️ No matching .gro for {itp_file}, skipping.")
        continue

    itp_atoms = get_atom_names_itp(itp_file)
    gro_atoms = get_atom_names_gro(gro_file)

    if itp_atoms != gro_atoms:
        mismatched.append(base)
        print(f"Atom name mismatch in {base}:")
        print(f"  ITP: {itp_atoms}")
        print(f"  GRO: {gro_atoms}")
    else:
        print(f"{base} atom names match.")

# Summary
print("\n--- Summary ---")
if mismatched:
    print("Mismatched files:", ", ".join(mismatched))
else:
    print("All .itp and .gro files have matching atom names.")

