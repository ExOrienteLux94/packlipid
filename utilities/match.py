import os
import re

def count_atoms_gro(filename):
    """Count atoms in a .gro file (second line gives atom count)."""
    with open(filename) as f:
        lines = f.readlines()
        if len(lines) < 3:
            return 0
        try:
            return int(lines[1].strip())
        except ValueError:
            return len(lines) - 3  # fallback if header line is non-numeric

def count_atoms_itp(filename):
    """Count atoms in a .itp file (count lines in [ atoms ] section)."""
    with open(filename) as f:
        lines = f.readlines()
    in_atoms = False
    count = 0
    for line in lines:
        line = line.strip()
        if line.startswith("["):
            in_atoms = line.startswith("[ atoms ]")
            in_atoms2 = line.startswith("[atoms]")
            continue
        if (in_atoms == True or in_atoms2 == True) and line and not line.startswith(";"):
            count += 1
    return count

# Compare matching files
gro_files = [f for f in os.listdir(".") if f.endswith(".gro")]
mismatched = []

for gro in gro_files:
    base = os.path.splitext(gro)[0]
    itp = base + ".itp"
    if os.path.exists(itp):
        gro_atoms = count_atoms_gro(gro)
        itp_atoms = count_atoms_itp(itp)
        if gro_atoms != itp_atoms:
            mismatched.append((base, gro_atoms, itp_atoms))

# Output results
if mismatched:
    print("Files with mismatched atom counts:")
    for base, g, i in mismatched:
        print(f"{base}: .gro={g}, .itp={i}")
else:
    print("All matching .gro and .itp files have consistent atom counts.")

