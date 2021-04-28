from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = automodel(env,
    alnfile  = 'prot.seq',
    knowns   = 'prot',
    sequence = 'prot_fill',
    assess_methods=(assess.DOPE,
                    assess.GA341))
a.starting_model= 1
a.ending_model  = 1
a.make()
a.write(file='out.pdb')

mdl = model(env, file='prot')
aln = alignment(env, file='prot.seq', align_codes=('prot', 'prot_fill'))
a.res_num_from(mdl,aln)
a.write(file='rename.pdb')
