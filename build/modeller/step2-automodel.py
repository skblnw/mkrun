from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(AutoModel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['F','H'],
                             renumber_residues=[1,1])

a = MyModel(env,
    alnfile  = 'prot.seq',
    knowns   = 'prot',
    sequence = 'prot_fill',
    assess_methods=(assess.DOPE,
                    assess.GA341))
a.starting_model= 1
a.ending_model  = 3

a.make()
