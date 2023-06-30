import pandas as pd
import numpy as np
import sys
import os
from scipy.stats import sem

# Read data from the file and create a DataFrame
def read_data(file_path):
    data_list = []
    with open(file_path, "r") as file:
        for line in file:
            elements = line.strip().split()
            if len(elements) == 5: 
                data_list.append(float(elements[1]))
            else:
                print(f'Error: Incomplete line encountered in file {file_path} This line is skipped.')
    data_array = np.array(data_list)
    return data_array.mean(), sem(data_array), len(data_array) # Return mean, standard error and count

def mutate_ref(ref, pos):
    if ref[pos - 1] != 'A':
        return ref[:pos-1] + 'A' + ref[pos:]
    else:
        print(f'Error: Original amino acid at position {pos} is already A')
        return None

def main(original_seq):
    # Define positions and states
    positions = [f"pos{i}" for i in range(1, len(original_seq)+1)]
    states = ['bound', 'free']

    # Initialize the list to store results
    results = []

    for position in positions:
        row = {}
        target_seq = mutate_ref(original_seq, int(position[-1])) # Mutate the reference sequence
        if target_seq is None:
            continue

        row['original'] = original_seq
        row['target'] = target_seq

        for state in states:
            file_path = os.path.join(position, state, "decompose_summary.dat")
            mean, error, count = read_data(file_path)
            row[f'{state}_dg'] = round(mean, 2)
            row[f'{state}_error'] = round(error, 2)
            row[f'{state}_count'] = count
        
        # Calculate ddG and its error
        row['ddG'] = round(row['bound_dg'] - row['free_dg'], 2)
        row['ddG_error'] = round(np.sqrt(row['bound_error']**2 + row['free_error']**2), 2)

        results.append(row)

    # Convert results to DataFrame and save to CSV
    df_results = pd.DataFrame(results)
    df_results.to_csv('results.csv', index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide the original sequence as an argument (e.g., 'IMDQVPFSV')")
    else:
        original_seq = sys.argv[1]
        main(original_seq)
