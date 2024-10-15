# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:51:25 2022

if or not strcutures can be predicted by NN 
y is the energy 
@author: jlchen
"""

import math
import numpy as np
import os
#from ase.calculators.reann import REANN
from ase.calculators.eann import EANN
#from ase.calculators.manual_written import Manual_written
from ase.io.trajectory import Trajectory
from ase.db import connect
from scipy import signal

import multiprocessing as mp
from concurrent import futures
import time
import scipy.signal as signal


class NN_predicted:

    def __init__(self, filename, filetype, atomtype1, atomtype2,
                 nn1='EANN_PES_DOUBLE1.pt', nn2='EANN_PES_DOUBLE1.pt'):
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype1
        self.nn1 = nn1
        self.nn2 = nn2
        self.filename = filename
        self.filetype = filetype
        self.calc1 = EANN(atomtype=atomtype1, nn=nn1)
        self.calc2 = EANN(atomtype=atomtype2, nn=nn2)

    def structure_energy(self):
        if self.filetype == 'traj':
            da = Trajectory(self.filename)
            structure = [da[i] for i in range(len(da))]
        elif self.filetype == 'db':
            da = connect(self.filename)
            structure = [da.get(id=i).toatoms() for i in range(1, len(da)+1)]
        else:
            raise ValueError(
                'there is wrong filetype (db or traj based on ASE)')
        return structure

    def get_Energy(self):
        self.atom = self.calc1
        energy_all1 = []
        energy_all2 = []
        structures = self.structure_energy()
        for i in range(len(structures)):
            atoms = structures[i]
            atoms.calc = self.calc1
            energy1 = atoms.get_potential_energy()
            energy_all1.append(energy1)
            del atoms.calc
            atoms.calc = self.calc2
            energy2 = atoms.get_potential_energy()
            energy_all2.append(energy2)
        return energy_all1, energy_all2

    def get_Force(self):
        self.atom = self.calc1
        forces_all1 = []
        forces_all2 = []
        structures = self.structure_energy()
        for i in range(len(structures)):
            atoms = structures[i]
            atoms.calc = self.calc1
            forces1 = atoms.get_forces()
            forces_all1.append(forces1)
            del atoms.calc
            atoms.calc = self.calc2
            forces2 = atoms.get_forces()
            forces_all2.append(forces2)
            print("#--------%s is finished-------------" % i)
        return forces_all1, forces_all2

    def energy_convergence(self, energy_min=0.02, energy_max=2):
        energy_all1, energy_all2 = self.get_Energy()
        structure_train = []
        structure_data = self.structure_energy()
        for i in range(len(energy_all1)):
            energy1 = energy_all1[i]
            energy2 = energy_all2[i]
            if abs(energy1-energy2) >= energy_min \
                    and abs(energy1-energy2) <= energy_max:
                structure_train.append(structure_dada[i])
        return structure_train

    def force_convergence(self, force_min=0.2, force_max=2):
        forces_all1, forces_all2 = self.get_Force()
        structure_train = []
        structure_data = self.structure_energy()
        print(len(forces_all1))
        for i in range(len(forces_all1)):
            forces1 = forces_all1[i]
            forces2 = forces_all2[i]
            Force_F1 = np.array([f_e for f_e in forces1 if np.std(f_e) != 0])
            Force_F2 = np.array([f_e for f_e in forces2 if np.std(f_e) != 0])
            #std_forces = abs(np.std(Force_F1)-np.std(Force_F2))
            std_forces = abs(np.std(Force_F1-Force_F2))
            print(std_forces)
            if std_forces >= force_min \
                    and std_forces <= force_max:
                structure_train.append(structure_data[i])
        return structure_train

# -------------------------filename is db------------------
# -------------------------filetype is the db or traj(default)-----------


class NSDS_algorithm:
    def __init__(self, filename,  atomtype1, atomtype2,
                 nn1='EANN_PES_DOUBLE.pt', nn2='EANN_PES_DOUBLE2.pt',
                 filetype='traj',
                 period1=[1, 1, 1], period2=[1, 1, 1], beta=2,
                 energy_criterion=2, distance=2):
        self.filename = filename
        self.filetype = filetype
        self.atomtype1 = atomtype1
        self.atomtype2 = atomtype1
        self.nn1 = nn1
        self.nn2 = nn2
        self.period1 = period1
        self.period2 = period2
        self.beta = beta
        if nn1 != None:
            self.calc1 = EANN(atomtype=atomtype1, period=period1, nn=self.nn1)
        self.calc2 = EANN(atomtype=atomtype2, period=period2, nn=self.nn2)
        self.energy_criterion = energy_criterion
        self.distance = distance

    # -----------add energy for NN2---------------
    def structure_energy(self):
        if self.filetype == 'traj':
            da = Trajectory(self.filename)
            energy = [da[i].get_potential_energy() for i in range(len(da))]
            structure = [da[i] for i in range(len(da))]
        elif self.filetype == 'db':
            da = connect(self.filename)
            energy = [da.get(id=i).energy for i in range(1, len(da)+1)]
            structure = [da.get(id=i).toatoms() for i in range(1, len(da)+1)]
        else:
            raise ValueError(
                'there is wrong filetype (db or traj based on ASE)')
        return structure, energy

    # ----------------use reann2-----------
    def calc_nn(self, atoms):
        atoms.calc = self.calc2
        en = atoms.get_potential_energy()
        print(en)
        return en

    def energy2_calc(self, structure):
        if not os.path.exists('./energy2.npy'):
            #    energy2 = []
            #    with futures.ThreadPoolExecutor(max_workers=processes) as executor:  # thread
            #     with futures.ProcessPoolExecutor(max_workers=processes) as executor: #process
            #       for en in executor.map(self.calc_nn,structure):
            #           energy2.append(en)
            energy2 = [self.calc_nn(atoms) for atoms in structure]
            np.save('energy2.npy', np.array(energy2))
        else:
            energy2 = np.load('energy2.npy')
        return energy2

    # #---------F = -w(y1-y2)-----------
    def NSDS(self, structure, energy1):
        # structure, energy1 = self.structure_energy()
        processes = min(mp.cpu_count(), len(structure))
        #print([s.get_positions()[-1] for s in structure])
        energy2 = self.energy2_calc(structure)
        # if not os.path.exists('./energy2.npy'):
        #    energy2 = []
        #    print(len(energy1),'\n',len([s.symbols for s in structure]))
        #    with futures.ThreadPoolExecutor(max_workers=processes) as executor:  # thread
        #     with futures.ProcessPoolExecutor(max_workers=processes) as executor: #process
        #       for en in executor.map(self.calc_nn,structure):
        #           energy2.append(en)
        #    energy2 = [self.calc_nn(atoms) for atoms in structure]
        #    np.save('energy2.npy',np.array(energy2))
        # else:
        #    energy2 = np.load('energy2.npy')
        energy1_m = np.array(energy1) - min(energy1)
        energy2_m = np.array(energy2) - min(energy1)
        for i in range(100):
            if self.updata_beta() == True:
                break
            else:
                if i >= 99:
                    raise ValueError('there is not appropriate beta')
        w = self.w_solve(energy1_m, energy2_m)
        F = -w*((np.array(energy1_m)-np.array(energy2_m))**2)
        return F

    # #--------w = e(-beta(y1**2+y2**2)/2)------------------
    def w_solve(self, y1, y2):
        y_mean = (np.array(y1)**2+np.array(y2)**2)/2
        w = math.e**(-self.beta*y_mean)
        return w

    # ---------------calc point (polyfit, we set 100)-------------------
    def extreme_point(self):
        structure, energy1 = self.structure_energy()
        energy1_m = np.array(energy1) - min(energy1)
        F = self.NSDS(structure, energy1)
        variable_coefficient = np.std(F)/np.mean(np.abs(F))
        np.save('F.npy', F)
        # print(variable_coefficient,np.std(F),np.mean(np.abs(F)))
        if abs(variable_coefficient) <= 0.1:
            raise ValueError(
                '-----------there is not data to train in %s energy criterion----------------' % self.energy_criterion)
        else:
            max_poly = min(len(F)/10, 100)
            # ------------------poly------------
            p1 = np.polyfit(energy1_m, F, max_poly)
            yvals = np.polyval(p1, energy1_m)
            peak = signal.find_peaks(yvals, distance=self.distance)[0]
            #prominences = signal.peak_prominences(yvals, peaks)[0]
            index = list(peak)
        # ----------------distance--------------------
        #index = signal.argrelextrema(F, np.greater, order=min(int(len(F)/10),5))
            energy2 = self.energy2_calc(structure)
            energy1_train = [energy1[i] for i in index]
            energy2_train = [energy2[i] for i in index]
            for i in range(len(energy1_train)):
                if abs(energy1_train[i]-energy2_train[i]) <= 0.02:
                    # print(index,energy1.index(energy1_train[i]))
                    # print(abs(energy1_train[i]-energy2_train[i]))
                    index.remove(energy1.index(energy1_train[i]))
            if len(index) != 0:
                structure_will_train = [structure[i] for i in index]
            else:
                raise ValueError(
                    '-----------there is not data to train in %s energy criterion----------------' % self.energy_criterion)

        return structure_will_train

    # ---------update beta for control the meta-energy-----------
    def updata_beta(self):
        w = math.e**(-self.beta*self.energy_criterion)
        if abs(w) <= 1E-2:
            return True
        else:
            self.beta = 1.1*self.beta
            return False
