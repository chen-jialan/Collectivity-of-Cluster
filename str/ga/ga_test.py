from random import random, randint
#import ramdom
from ase.io import write
from ase.optimize import LBFGS
#from ase.calculators.emt import EMT
from ase.calculators.eann import EANN
from ase.calculators.manual_written import Manual_written
import os
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


def calc(a):
    atomtype = ['Cu', 'O', 'C']
    a.calc = EANN(atomtype=atomtype, period=[1, 1, 1], nn='EANN_PES_DOUBLE.pt')
    # dyn = LBFGS(a, trajectory=None, logfile=None)
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


if __name__ == "__main__":
    # Change the following three parameters to suit your needs
    population_size = 48
    mutation_probability = 0.5
    n_to_test = 200
    processes = mp.cpu_count()
    if os.path.exists('./tmp_traj'):
        os.system('rm -r ./tmp_traj')

    # Initialize the different components of the GA
    da = DataConnection('gadb.db')
    atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
    n_to_optimize = len(atom_numbers_to_optimize)
    slab = da.get_slab()
    all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
    blmin = closest_distances_generator(all_atom_types,
                                        ratio_of_covalent_radii=0.7)

    comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                         pair_cor_cum_diff=0.015,
                                         pair_cor_max=0.7,
                                         dE=0.02,
                                         mic=False)
    # comp = OFPComparator(n_top=n_to_optimize,
    #                     cos_dist_max=0.01)
    pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
    mutations = OperationSelector([1., 15., 1., 5.],
                                  [MirrorMutation(blmin, n_to_optimize),
                                   RattleMutation(
                                       blmin, n_to_optimize, rattle_strength=2.3),
                                   PermutationMutation(n_to_optimize),
                                   RotationalMutation(blmin, n_to_optimize)])

    # Relax all unrelaxed structures (e.g. the starting population)
    number_da = da.get_number_of_unrelaxed_candidates()
    if number_da > 0:
        print('Relaxing starting candidate')
        a = da.get_all_unrelaxed_candidates()

        with futures.ThreadPoolExecutor(max_workers=number_da) as executor:  # thread
            # with futures.ProcessPoolExecutor(max_workers=processes) as executor: #process
            for a_relax in executor.map(calc, a):
                da.add_relaxed_step(a_relax)
        print('Relaxing done candidate')
    # create the population
    population = Population(data_connection=da,
                            population_size=population_size,
                            comparator=comp, use_extinct=True)
    cc = GenerationRepetitionConvergence(population, number_of_generations=3,
                                         number_of_individuals=-1, max_generations=n_to_test)

    # # test n_to_test new candidates
    for i in range(n_to_test):
        if cc.converged():
            print('converged')
            break
        print('Now starting configuration number {0}'.format(i))
        offspring = []
        # for i in range(randint(1,processes)):
        for j in range(processes):
            a1, a2 = population.get_two_candidates()
            a3, desc = pairing.get_new_individual([a1, a2])
            if a3 is None:
                continue
            da.add_unrelaxed_candidate(a3, description=desc)

            # Check if we want to do a mutation
            if random() < mutation_probability:
                a3_mut, desc = mutations.get_new_individual([a3])
                if a3_mut is not None:
                    da.add_unrelaxed_step(a3_mut, desc)
                    a3 = a3_mut
            offspring.append(a3)
        print('Now starting offsprings %s' % len(offspring))
        # thread
        with futures.ThreadPoolExecutor(max_workers=len(offspring)) as executor:
            # with futures.ProcessPoolExecutor(max_workers=processes) as executor: #process
            for a_relax in executor.map(calc, offspring):
                da.add_relaxed_step(a_relax)
                population.update()

        if i % 10 == 9:
            if not os.path.exists('./tmp_traj'):
                os.mkdir('./tmp_traj')
            write('all_candidates.traj', da.get_all_relaxed_candidates())
            os.system(
                'mv all_candidates.traj tmp_traj/all_candidates_%s.traj' % i)

    write('all_candidates.traj', da.get_all_relaxed_candidates())
