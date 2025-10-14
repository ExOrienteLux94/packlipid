import os
import re
from collections import defaultdict

# Regex to find molecule names after [ moleculetype ]
molecule_pattern = re.compile(r'^\s*\[\s*moleculetype\s*\]\s*\n\s*([^\s;]+)', re.M)

file_molecules = {}
all_molecules = defaultdict(list)
intra_duplicates = []

# Process all .itp files
for fname in sorted(f for f in os.listdir('.') if f.endswith('.itp')):
    with open(fname) as f:
        lines = f.readlines()

    # Remove comment lines
    lines = [line for line in lines if not line.strip().startswith(';')]

    text = "".join(lines)
    molecules = molecule_pattern.findall(text)

    # Check duplicates within a file
    seen = set()
    dups = set()
    for mol in molecules:
        if mol in seen:
            dups.add(mol)
        seen.add(mol)

    if dups:
        intra_duplicates.append((fname, sorted(dups)))

    file_molecules[fname] = molecules

    # Track globally for cross-file duplicates
    for mol in molecules:
        all_molecules[mol].append(fname)

# Cross-file duplicates
cross_duplicates = {mol: files for mol, files in all_molecules.items() if len(files) > 1}

# --- Reporting ---
print("=== Molecule Name Consistency Check ===\n")

if intra_duplicates:
    print("❌ Duplicate molecule names *within* files:")
    for fname, dups in intra_duplicates:
        print(f"  {fname}: {', '.join(dups)}")
else:
    print("✅ No intra-file duplicates found.")

print("\n---")

if cross_duplicates:
    print("❌ Duplicate molecule names *across* files:")
    for mol, files in sorted(cross_duplicates.items()):
        print(f"  Molecule '{mol}' found in: {', '.join(sorted(files))}")
else:
    print("✅ No cross-file duplicates found.")

print("\n--- Done ---")

