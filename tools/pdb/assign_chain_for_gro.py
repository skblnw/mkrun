import sys
import argparse

def assign_chain_ids(filename, chains, starts):
    with open(filename, 'r') as f:
        lines = f.readlines()

    prev_res_id = None
    chain_counter = 0
    new_lines = []

    for line in lines:
        if line.startswith('ATOM'):
            res_id = int(line[22:26].strip())

            if prev_res_id is not None and (res_id != prev_res_id and res_id != prev_res_id + 1):
                chain_counter += 1

            chain_id = chains[chain_counter] if chain_counter < len(chains) else 'X'
            resid_offset = starts[chain_counter]

            # Adjust the residue ID by subtracting the initial residue id of the previous chain and adding the starting residue id of the current chain
            new_res_id = res_id + resid_offset - 1

            # Update the line with the new chain id and the new residue ID
            new_line = line[:21] + chain_id + f"{new_res_id:>4}" + line[26:]
            new_lines.append(new_line)

            prev_res_id = res_id

        else:
            new_lines.append(line)

    with open('new_' + filename, 'w') as f:
        f.writelines(new_lines)

def main():
    parser = argparse.ArgumentParser(description="Assign specified chain IDs and start residue numbers to a PDB file.")
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("--chains", type=str, nargs='+', help="List of chain IDs to assign.", required=True)
    parser.add_argument("--start", type=int, nargs='+', help="List of starting residue numbers for each chain.", required=True)

    args = parser.parse_args()

    if len(args.chains) != len(args.start):
        print("Error: The number of chains should match the number of starting residue numbers.")
        exit(1)

    assign_chain_ids(args.pdb_file, args.chains, args.start)

if __name__ == "__main__":
    main()
