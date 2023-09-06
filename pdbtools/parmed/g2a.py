#!/usr/bin/env python3

import parmed as pmd
import argparse


def do_pmd(topology, structure, coord_only):

    gmx = pmd.load_file(topology, xyz=structure)
    if not coord_only:
        gmx.save('ionized.prmtop')
    gmx.save('out.ncrst')

def main(args):

    do_pmd(args.topology, args.structure, args.coord_only)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processing PDB and PSF files.")
    parser.add_argument("topology", type=str, help="Path to the Topology file.")
    parser.add_argument("structure", type=str, help="Path to the PDB or GRO file.")
    parser.add_argument("--coord-only", action='store_true', help="Save coordinates only.")
    
    args = parser.parse_args()
    main(args)
