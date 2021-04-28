from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

class MyModel(automodel):
    num_map = None

    def fix_numbering(self):
        for old_resi, new_resi in zip(self.residues, self.num_map):
            old_resi.num = str(new_resi)

    def user_after_single_model(self):
        self.fix_numbering()

env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = MyModel(env,
    alnfile  = 'prot.seq',
    knowns   = 'prot',
    sequence = 'prot_fill',
    assess_methods=(assess.DOPE,
                    assess.GA341))
a.num_map = [..., 1400, 1401, 1450, 1451,...]

a.starting_model= 1
a.ending_model  = 1
a.make()
a.write(file='out.pdb')

mdl = model(env, file='prot')
aln = alignment(env, file='prot.seq', align_codes=('prot', 'prot_fill'))
a.res_num_from(mdl,aln)
a.write(file='rename.pdb')
