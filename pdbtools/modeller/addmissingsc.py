import sys
from modeller import *
from modeller.automodel import *

# Check if a PDB file was provided as an argument
if len(sys.argv) != 2:
    print("Usage: python addmissingsc.py pdb_file.pdb")
    sys.exit(1)

pdb_file = sys.argv[1]
pdb_code = pdb_file[:-4]
sequence_file = f"{pdb_code}.seq"
alignment_file = f"{pdb_code}.ali"

# Step 1: Extract the sequence from the PDB file
log.verbose()
env = environ()
env.io.atom_files_directory = ['.', '/path/to/your/pdb/files']
mdl = model(env, file=pdb_file)
aln = alignment(env)
aln.append_model(mdl, align_codes=pdb_file)
aln.write(file=sequence_file)

# Step 2: Create an alignment file by duplicating the sequence
with open(sequence_file) as seq_file:
    lines = seq_file.readlines()
    sequence = ''.join(lines[2:])

with open(alignment_file, 'w') as ali_file:
    ali_file.write(f">P1;{pdb_code}\n")
    ali_file.write(sequence)
    ali_file.write(f">P1;{pdb_code}_fill\n")
    ali_file.write(sequence)

# Step 3: Load the structure, identify missing side chains, and add them
class AddMissingSideChains(automodel):
    def select_atoms(self):
        # Here, we're returning all atoms, assuming we want to check the entire structure
        # for missing side chains. Adjust if needed.
        return selection(self)

# Load the alignment file and known structure
a = alignment(env, file=alignment_file)
m = AddMissingSideChains(env, alnfile=alignment_file,
                         knowns=pdb_code, sequence=pdb_code, assess_methods=(assess.DOPE,assess.GA341))
m.starting_model = 1
m.ending_model = 1

# Perform the task
m.very_fast()
m.make()
