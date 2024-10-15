from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.constraints import FixAtoms
import numpy as np
from ase.build import fcc111
import os
import ase.io.vasp

db_file = 'gadb.db'
if os.path.exists('./%s' % db_file):
    os.remove('./%s' % db_file)
# create the surface
slab = ase.io.vasp.read_vasp("POSCAR_base")
pos = slab.get_positions()
cell = slab.get_cell()
p0 = np.array([min(pos[:, 0]), min(pos[:, 1])/3, max(pos[:, 2])+1.6])
v1 = cell[0, :] * 1
v2 = cell[1, :] * 1
v3 = cell[2, :]
print(v1, v2, v3)
v3[2] = 2

# Define the composition of the atoms to optimize,47 is Ag, 79 is Au
atom_numbers = 10 * [8] + 8 * [29]

# define the closest distance two atoms of a given species can be to each other
unique_atom_types = get_all_atom_types(slab, atom_numbers)
blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=0.75)

# create the starting population
sg = StartGenerator(slab, atom_numbers, blmin,
                    box_to_place_in=[p0, [v1, v2, v3]])

# generate the starting population
population_size = 24
starting_population = [sg.get_new_candidate() for i in range(population_size)]

# create the database to store information in
d = PrepareDB(db_file_name=db_file,
              simulation_cell=slab,
              stoichiometry=atom_numbers)

for a in starting_population:
    d.add_unrelaxed_candidate(a)
