from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

class MyModel(loopmodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A', 'B'],
                             renumber_residues=[1, 1])
        # # Another way to label individual chains:
        # self.chains[0].name = 'A'
        # self.chains[1].name = 'B'
		# self.chains[2].name = 'C'

env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = MyModel(env,
				alnfile  = 'alignment.ali',
				knowns   = 'chainA',
				sequence = 'chainA_loop')
a.starting_model= 1
a.ending_model  = 1

a.loop.starting_model = 1
a.loop.ending_model   = 2
a.loop.md_level       = refine.fast

a.make()