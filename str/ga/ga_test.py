import os
import argparse
from random import random, randint
from ase.io import write
from ase.optimize import LBFGS
from ase.calculators.eann import EANN
from ase.calculators.manual_written import Manual_written
from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
from ase.ga.standardmutations import RotationalMutation
from ase.ga.convergence import GenerationRepetitionConvergence
from ase.ga.ofp_comparator import OFPComparator
import multiprocessing as mp
from concurrent import futures


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run genetic algorithm for atomic structure optimization.")
    parser.add_argument("--population_size", type=int,
                        default=48, help="Size of the population.")
    parser.add_argument("--mutation_probability", type=float,
                        default=0.5, help="Probability of mutation.")
    parser.add_argument("--n_to_test", type=int, default=200,
                        help="Number of configurations to test.")
    parser.add_argument("--processes", type=int, default=os.cpu_count(),
                        help="Number of processes or threads.")
    parser.add_argument("--db_file", type=str,
                        default="gadb.db", help="Database file for ASE-GA.")
    parser.add_argument("--atom_types", nargs="+", default=[
                        "Cu", "O", "C"], help="Atomic types for the EANN calculator (e.g., Cu O C)")
    parser.add_argument("--cutoff_radius", type=float,
                        default=0.7, help="Covalent radii ratio cutoff.")
    return parser.parse_args()


def calc(a):
    """Relax atoms and compute energy/forces."""
    atomtype = ["Cu", "O", "C"]
    a.calc = EANN(atomtype=atomtype, period=[1, 1, 1], nn='EANN_PES_DOUBLE.pt')
    dyn = LBFGS(a)
    dyn.run(fmax=0.1, steps=200)
    energy = a.get_potential_energy()
    force = a.get_forces()
    a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()
    del a.calc
    a.calc = Manual_written(energy_write=energy, force_write=force)
    e = a.get_potential_energy()
    f = a.get_forces()
    return a


def setup_ga(args):
    """Set up the genetic algorithm components."""
    # Initialize the database and other components
    da = DataConnection(args.db_file)
    atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
    n_to_optimize = len(atom_numbers_to_optimize)
    slab = da.get_slab()

    # Setup the unique atom types and bond length parameters
    all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
    blmin = closest_distances_generator(
        all_atom_types, ratio_of_covalent_radii=args.cutoff_radius)

    comp = InteratomicDistanceComparator(
        n_top=n_to_optimize, pair_cor_cum_diff=0.015, pair_cor_max=0.7, dE=0.02, mic=False)
    pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)

    mutations = OperationSelector(
        [1., 15., 1., 5.],
        [
            MirrorMutation(blmin, n_to_optimize),
            RattleMutation(blmin, n_to_optimize, rattle_strength=2.3),
            PermutationMutation(n_to_optimize),
            RotationalMutation(blmin, n_to_optimize)
        ]
    )

    # Setup population
    population = Population(
        data_connection=da, population_size=args.population_size, comparator=comp, use_extinct=True)
    cc = GenerationRepetitionConvergence(
        population, number_of_generations=3, number_of_individuals=-1, max_generations=n_to_optimize)

    return da, population, cc, pairing, mutations


def relax_starting_population(da, processes):
    """Relax starting population candidates."""
    number_da = da.get_number_of_unrelaxed_candidates()
    if number_da > 0:
        print('Relaxing starting candidate')
        candidates = da.get_all_unrelaxed_candidates()

        with futures.ThreadPoolExecutor(max_workers=number_da) as executor:
            for a_relax in executor.map(calc, candidates):
                da.add_relaxed_step(a_relax)
        print('Relaxing done candidate')


def generate_offspring(population, pairing, mutations, mutation_probability, processes):
    """Generate offspring by crossover and mutation."""
    offspring = []
    for j in range(processes):
        a1, a2 = population.get_two_candidates()
        a3, desc = pairing.get_new_individual([a1, a2])
        if a3 is None:
            continue
        offspring.append(a3)
        if random() < mutation_probability:
            a3_mut, desc = mutations.get_new_individual([a3])
            if a3_mut is not None:
                offspring.append(a3_mut)
    return offspring


def main():
    args = parse_args()

    # Setup GA components
    da, population, cc, pairing, mutations = setup_ga(args)

    # Relax the starting population
    relax_starting_population(da, args.processes)

    # GA iterations
    for i in range(args.n_to_test):
        if cc.converged():
            print('Converged!')
            break
        print(f"Now starting configuration number {i}")

        # Generate offspring
        offspring = generate_offspring(
            population, pairing, mutations, args.mutation_probability, args.processes)

        # Relax offspring and update population
        with futures.ThreadPoolExecutor(max_workers=len(offspring)) as executor:
            for a_relax in executor.map(calc, offspring):
                da.add_relaxed_step(a_relax)
                population.update()

        # Output and save candidates
        if i % 10 == 9:
            write('all_candidates.traj', da.get_all_relaxed_candidates())
            os.system(
                f'mv all_candidates.traj tmp_traj/all_candidates_{i}.traj')

    write('all_candidates.traj', da.get_all_relaxed_candidates())


if __name__ == "__main__":
    main()
