from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(AutoModel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['A'],
                             renumber_residues=[69])
    def select_atoms(self):
        return Selection(self.residue_range('400:A', '449:A'),
                         self.residue_range('822:A', '836:A'))

a = MyModel(env,
    alnfile  = 'prot.seq',
    knowns   = 'prot',
    sequence = 'prot_fill',
    assess_methods=(assess.DOPE,
                    assess.GA341))
a.starting_model= 1
a.ending_model  = 3

a.make()
