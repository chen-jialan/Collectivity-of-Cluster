from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
import numpy as np
import os
import ase.io.vasp

import os
import argparse
import numpy as np
from ase import io
from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator, get_all_atom_types


def expand_atom_spec(atom_specs):
    """
    Convert list like ['8x10', '29x8'] to [8, 8, ..., 29, 29]
    """
    expanded = []
    for spec in atom_specs:
        if 'x' in spec:
            z, n = spec.split('x')
            expanded.extend([int(z)] * int(n))
        else:
            expanded.append(int(spec))
    return expanded


def generate_population(
    poscar_file,
    db_file,
    atom_numbers,
    population_size,
    covalent_ratio=0.75,
    vacuum_height=1.6
):
    # Load slab structure
    slab = io.read(poscar_file)
    pos = slab.get_positions()
    cell = slab.get_cell()

    # Define region to place atoms
    p0 = np.array([min(pos[:, 0]), min(pos[:, 1]), max(pos[:, 2]) + vacuum_height])
    v1, v2, v3 = cell[0], cell[1], [0, 0, 2.0]  # Thin z region

    # Generate minimum bond length dictionary
    unique_atom_types = get_all_atom_types(slab, atom_numbers)
    blmin = closest_distances_generator(
        atom_numbers=unique_atom_types,
        ratio_of_covalent_radii=covalent_ratio
    )

    # Start generator
    sg = StartGenerator(
        slab,
        atom_numbers,
        blmin,
        box_to_place_in=[p0, [v1, v2, v3]]
    )

    # Generate population
    population = [sg.get_new_candidate() for _ in range(population_size)]

    # Prepare database
    if os.path.exists(db_file):
        os.remove(db_file)

    db = PrepareDB(
        db_file_name=db_file,
        simulation_cell=slab,
        stoichiometry=atom_numbers
    )

    for a in population:
        db.add_unrelaxed_candidate(a)

    print(f"✅ Generated {population_size} candidates into {db_file}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate GA initial population.")
    parser.add_argument("--poscar", type=str, default="POSCAR_base", help="Input POSCAR file")
    parser.add_argument("--db", type=str, default="gadb.db", help="Output .db file for ASE-GA")
    parser.add_argument("--atoms", nargs="+", default=["8x10", "29x8"],
                        help="Atoms and counts, e.g., 8x10 29x8 for 10 O and 8 Cu")
    parser.add_argument("--size", type=int, default=24, help="Population size")
    parser.add_argument("--zvac", type=float, default=1.6, help="Vacuum above surface for placement")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    atom_numbers = expand_atom_spec(args.atoms)

    generate_population(
        poscar_file=args.poscar,
        db_file=args.db,
        atom_numbers=atom_numbers,
        population_size=args.size,
        vacuum_height=args.zvac
    )
