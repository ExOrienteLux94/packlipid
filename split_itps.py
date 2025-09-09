#!/usr/bin/env python3
import re
import sys
from pathlib import Path

def parse_defines(ffbonded_text):
    """Parse #define macros into dictionary."""
    define_pattern = re.compile(r"^[ \t]*#define[ \t]+(\S+)[ \t]+(.+?)\s*$", re.MULTILINE)
    return {m.group(1): m.group(2) for m in define_pattern.finditer(ffbonded_text)}

def expand_macros(text, defines):
    """Replace macro tokens with their values using longest-name-first strategy."""
    keys_sorted = sorted(defines.keys(), key=len, reverse=True)
    for key in keys_sorted:
        pattern = re.compile(rf"\b{re.escape(key)}\b")
        text = pattern.sub(defines[key], text)
    return text

def split_molecules(expanded_text, output_dir):
    """Split expanded itp into separate files by moleculetype name."""
    lines = expanded_text.splitlines()
    blocks = []
    current_block = []
    current_name = None
    in_molecule = False
    save_lines = False

    for line in lines:
        # Start of new moleculetype
        #print(line)
        if line.strip().lower().startswith("[moleculetype") or line.strip().lower().startswith("[ moleculetype"):
            #print('lol')
            # Save the previous block
            if in_molecule and current_name and current_block:
                blocks.append((current_name, "\n".join(current_block)))
            # Reset
            current_block = [line]
            current_name = None
            in_molecule = True
            save_lines = True
        elif save_lines:
            #print("yay")
            current_block.append(line)
            # Capture molecule name: first non-comment, non-empty line after [ moleculetype ]
            if in_molecule and current_name is None:
                if line.strip() and not line.strip().startswith(";") and not line.strip().startswith("["):
                    current_name = line.split()[0]
                    print(current_name)

    # Save the last block
    if in_molecule and current_name and current_block:
        blocks.append((current_name, "\n".join(current_block)))

    # Write each molecule to a file
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, block in blocks:
        outpath = output_dir / f"{name}.itp"
        with open(outpath, "w") as f:
            f.write(block + "\n")
        print(f"Written: {outpath}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python expand_split_itp.py martini_v3.0.0_ffbonded_v2.itp martini_v3.0.0_phospholipids_PC_v2.itp [output_dir]")
        sys.exit(1)

    ffbonded_file = sys.argv[1]
    lipid_itp_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "expanded_split"

    # Read files
    with open(ffbonded_file, "r") as f:
        ffbonded_text = f.read()
    with open(lipid_itp_file, "r") as f:
        lipid_text = f.read()

    # Expand macros
    defines = parse_defines(ffbonded_text)
    expanded_text = expand_macros(lipid_text, defines)

    # Write combined expanded file
    combined_path = Path(output_dir) / (Path(lipid_itp_file).stem + "_expanded.itp")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    with open(combined_path, "w") as f:
        f.write(expanded_text)
    print(f"Written combined expanded file: {combined_path}")

    # Split into individual files
    split_molecules(expanded_text, output_dir)

if __name__ == "__main__":
    main()
