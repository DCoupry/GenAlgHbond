from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
import numpy as np
from ase.build import molecule
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read 

db_file = 'gadb.db'

# create the surface
mof  = read("mil53.cif")
mof.set_constraint(FixAtoms(mask=len(mof) * [True]))
# define the volume in which the adsorbed water is optimized
# the volume is defined by a corner position (p0)
# and three spanning vectors (v1, v2, v3)
pos   = np.zeros(3)
cell  = mof.get_cell()

# Define the composition of the atoms to optimize
atom_numbers = 6*[29]

# define the closest distance two atoms of a given species can be to each other
unique_atom_types = get_all_atom_types(mof, atom_numbers)
cd = closest_distances_generator(atom_numbers=unique_atom_types,
                                 ratio_of_covalent_radii=1.2)

# create the starting population
sg = StartGenerator(slab=mof,
                    atom_numbers=atom_numbers,
                    closest_allowed_distances=cd,
                    box_to_place_in=[pos, cell])

# generate the starting population
population_size = 20
starting_population = []
for i in range(population_size):
	print ("POP+1")
	starting_population.append(sg.get_new_candidate() )

from ase.visualize import view   # uncomment these lines
view(starting_population)        # to see the starting population

# create the database to store information in
d = PrepareDB(db_file_name=db_file,
              simulation_cell=mof,
              stoichiometry=atom_numbers)

for a in starting_population:
    d.add_unrelaxed_candidate(a)
