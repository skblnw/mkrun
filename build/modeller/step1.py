from modeller import *
# Get the sequence of the 1qg8 PDB file, and write to an alignment file
code = 'prot'

env = Environ()
mdl = Model(env, file=code)
aln = Alignment(env)
aln.append_model(mdl, align_codes=code)
aln.write(file=code+'.seq')
