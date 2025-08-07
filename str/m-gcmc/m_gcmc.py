# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 10:13:34 2021
T = 350K, pressure = 1P
O:20(base on O2/2) CO:20
radio_CO = 0.5
CO=-.14802542E+02 O2=-.98551638E+01
350K:CO=-0.806 O2=-0.723 eV
300K:CO=-0.659 O2=-0.588 eV
#150K:CO=-0.233 O2=-0.197 eV
200K:CO=-0.417 O2=-0.432 eV 
100K:CO=-0.232 O2=-0.240 eV
400K:CO=-0.827 O2=-0.858 eV
@author: DELL
"""

# import math
from math import sqrt, exp, pi, cos, sin, log
import os
import shutil
import re
import numpy as np
from random import uniform, shuffle, sample, randint, choice
import linecache
import time
import copy
from devide_arc import mkdir_file
from POSCAR_to_atom import POSCAR_To_atom1, Energy_VASP
import copy
from tranfer import write_POSCAR, Energy_lammps, lammps2atom, badstr_energy
#from BFGS_method import bfgs_main
from ase.optimize import LBFGS
from ase.calculators.eann import EANN
from ase import Atoms
from ase.io.trajectory import Trajectory
import ase.io.vasp
# import dire2cart


class Population:
    def __init__(self, ratio_CO=0.5, adjust_value=0, pressure=1e5, 
                 standard_pressure=100000, ads_number=10, cycle_number=1, 
                 ConstantR=8.618E-5, temperature=300, cutoff_distance=1.95, 
                 cutoff_distance1=1.2, cutoff_distance2=2.7, cutoff_distance3=1.3, 
                 potential_O2=-1.071316E+01/2, potential_CO=-15.78612, 
                 potential_CO2=-24.24603, all_CO=20, all_O2=20):
        self.ratio_CO = ratio_CO
        self.ads_number = ads_number
        self.cycle_number = cycle_number
        self.ConstantR = ConstantR
        self.temperature = temperature
        self.cutoff_distance = cutoff_distance
        self.cutoff_distance1 = cutoff_distance1
        self.cutoff_distance2 = cutoff_distance2
        self.cutoff_distance3 = cutoff_distance3
        self.potential_O2 = potential_O2
        self.potential_CO = potential_CO
        self.potential_CO2 = potential_CO2
        self.all_CO = all_CO
        self.all_O2 = all_O2
        self.pressure = pressure
        self.standard_pressure = standard_pressure
        self.adjust_value = adjust_value

# -------------read all struture and position----------------------
    def atom(self, path, filename):
        list1 = list()
        file = path + '/' + filename
        with open(file, 'r') as fp:
            for line in fp:
                list1.append(line.strip())

        a = len(list1)
        at_dec = []
        atom = []
        for i in range(5, a-2):
            text = list1[i]
            line_all = re.findall(r"-?\d+\.?\d*e?-?\d*?", text)
            split1 = '%s-%s' % (str(text.split(' ', 1)[0]), i-4)
            atom.append(split1)

            x = float(line_all[0])
            y = float(line_all[1])
            z = float(line_all[2])

            at_dec.append([x, y, z])

        structure = dict(zip(atom, at_dec))
        return structure

# -------------read stuture on support----------------------
    def support_atom(self, structure, Height):
        atom_dis = structure
        dictO = dict()
        dictCe = dict()
        H = []
        for pot in atom_dis.keys():
            if 'Cu-' in pot:
                b = atom_dis[pot][2]
                H.append(b)
        b = min(H)

        for pot in atom_dis.keys():
            if 'O-' in pot:
                if float(atom_dis[pot][2]) <= (Height + 0.3):
                    dictO[pot] = atom_dis[pot]

        for pot in atom_dis.keys():
            if 'Ce-' in pot:
                dictCe[pot] = atom_dis[pot]

        support = dict()
        support.update(dictCe)
        support.update(dictO)
        height = []
        for pot in dictO.keys():
            height.append(dictO[pot][2])
        max_H = max(height)
        return support, max_H

# -------------read cluster on stuture----------------------
    def cluster_atom(self, structure, Height):
        atom_dis = structure
        dictCu = dict()
        dictO = dict()
        cluster_and_gas = dict()
        atom_dis_cluster = dict()
        cluster = dict()
        cluster_new = dict()
    # --------------confirm cluster&gas------------------
        for pot in atom_dis.keys():
            if 'Cu-' in pot:
                dictCu[pot] = atom_dis[pot]
                cluster[pot] = atom_dis[pot]
            else:
                if float(atom_dis[pot][2]) > Height:
                    cluster[pot] = atom_dis[pot]

        dis = []
        cluster_new.update(dictCu)

        cluster_and_gas_copy = copy.deepcopy(cluster_and_gas)
        for pot in cluster_and_gas.keys():
            if 'O-' in pot:
                arr_cluster = np.array(cluster_and_gas[pot])
                dis_min1, pot_min1 = self.min_distance(
                    arr_cluster, pot, cluster)
                cluster_1 = copy.deepcopy(cluster)
                del cluster_1[pot_min1]
                dis_min2, pot_min2 = self.min_distance(
                    arr_cluster, pot, cluster_1)
                if dis_min2 <= 2.3:
                    if 'O-' in pot_min1:
                        arr_cluster_post = np.array(cluster_and_gas[pot_min1])
                        dis_min1_post, pot_min1_post = self.min_distance(
                            arr_cluster_post, pot_min1, cluster)
                        cluster_1_post = copy.deepcopy(cluster)
                        del cluster_1_post[pot_min1_post]
                        dis_min2_post, pot_min2_post = self.min_distance(
                            arr_cluster_post, pot_min1, cluster_1_post)
                        if dis_min2_post <= 2.3:
                            cluster_new[pot] = cluster_and_gas[pot]
                    elif 'O-' in pot_min2:
                        arr_cluster_post = np.array(cluster_and_gas[pot_min2])
                        dis_min1_post, pot_min1_post = self.min_distance(
                            arr_cluster_post, pot_min2, cluster)
                        cluster_1_post = copy.deepcopy(cluster)
                        del cluster_1_post[pot_min1_post]
                        dis_min2_post, pot_min2_post = self.min_distance(
                            arr_cluster_post, pot_min2, cluster_1_post)
                        if dis_min2_post <= 2.3:
                            cluster_new[pot] = cluster_and_gas[pot]
                    elif 'Cu-' in pot_min1 and 'Cu-' in pot_min2:
                        cluster_new[pot] = cluster_and_gas[pot]

        return cluster_new

# -------------read stuture on gas----------------------
    def gas_atom(self, structure, support, cluster):
        atom_dis = copy.deepcopy(structure)
        support_dis = copy.deepcopy(support)
        cluster_dis = copy.deepcopy(cluster)
        for pot in support_dis.keys():
            del atom_dis[pot]
        for pot1 in cluster_dis.keys():
            del atom_dis[pot1]
        gas = atom_dis
        return gas

# --------------new cluster-move random----------------------------------------------------
    def new_cluster_move(self, cluster, H, allatom, cluster_gas_copy):
        allatom_dis = copy.deepcopy(allatom)
        cluster_dis = copy.deepcopy(cluster)
        cluster_gas = copy.deepcopy(cluster_gas_copy)
        disx = []
        disy = []
        disz = []
        for pot in cluster_dis.keys():
            disx.append(cluster_dis[pot][0])
            disy.append(cluster_dis[pot][1])
            disz.append(cluster_dis[pot][2])

        dlist = []
        number_min = len(cluster_dis)  # 8 is the number of metal cluster
        pot_all = []
        for i in range(1000):
            for pot_choice in cluster_dis.keys():
                if pot_choice not in dlist and uniform(0, 1) < 0.3:
                    dlist.append(pot_choice)
            if len(dlist) >= 1:
                break

        shuffle(dlist)
        for pot0 in dlist:
            if 'Cu' in pot0:
                arr_1 = cluster_dis[pot0]
                #arr_add = [uniform(-0.3,0.3),uniform(-0.3,0.3),uniform(-0.1,0.1)]
                arr_add = [random_2part(-0.3, -0.01, 0.01, 0.3), random_2part(-0.3, -
                                                                              0.01, 0.01, 0.3), random_2part(-0.3, -0.01, 0.01, 0.3)]
                arr1_1 = [arr_1[0]+arr_add[0], arr_1[1] +
                          arr_add[1], max(H, arr_1[2]+arr_add[2])]

                arr1 = np.array(arr1_1)

                move_dis = [0, 0, 0]
                arr1_1 = list(arr1)

                allatom_dis[pot0] = arr1_1
                cluster_dis[pot0] = arr1_1
                for pot_cluster in cluster_gas.keys():
                    if pot0 == pot_cluster:
                        if not 'none' in cluster_gas[pot_cluster]:
                            for pot_gas in cluster_gas[pot_cluster].keys():
                                x_gas = float(
                                    cluster_gas[pot_cluster][pot_gas][0])+arr_add[0]+move_dis[0]
                                y_gas = float(
                                    cluster_gas[pot_cluster][pot_gas][1])+arr_add[1]+move_dis[1]
                                z_gas = float(
                                    cluster_gas[pot_cluster][pot_gas][2])+arr_add[2]+move_dis[2]
                                allatom_dis[pot_gas] = [x_gas, y_gas, z_gas]
                                cluster_gas[pot_cluster][pot_gas] = [
                                    x_gas, y_gas, z_gas]
                # cluster_new,cluster_gas_new = self.single_atom(cluster_dis,cluster_gas)

            else:
                arr_1 = cluster_dis[pot0]
                #arr_add = [uniform(-0.3,0.3),uniform(-0.3,0.3),uniform(-0.1,0.1)]
                arr_add = [random_2part(-0.6, -0.1, 0.1, 0.6), random_2part(-0.6, -
                                                                            0.1, 0.1, 0.6), random_2part(-0.6, -0.1, 0.1, 0.6)]
                arr1_1 = [arr_1[0]+arr_add[0], arr_1[1] +
                          arr_add[1], max(H+0.1, arr_1[2]+arr_add[2])]

                arr1 = np.array(arr1_1)
                move_dis = [0, 0, 0]
                arr1_1 = list(arr1)

                allatom_dis[pot0] = arr1_1
                cluster_dis[pot0] = arr1_1
                for pot_cluster in cluster_gas.keys():
                    if pot0 == pot_cluster:
                        if not 'none' in cluster_gas[pot_cluster]:
                            for pot_gas in cluster_gas[pot_cluster].keys():
                                x_gas = float(
                                    cluster_gas[pot_cluster][pot_gas][0])+arr_add[0]+move_dis[0]
                                y_gas = float(
                                    cluster_gas[pot_cluster][pot_gas][1])+arr_add[1]+move_dis[1]
                                z_gas = float(
                                    cluster_gas[pot_cluster][pot_gas][2])+arr_add[2]+move_dis[2]
                                allatom_dis[pot_gas] = [x_gas, y_gas, z_gas]
                                cluster_gas[pot_cluster][pot_gas] = [
                                    x_gas, y_gas, z_gas]
#

# -------------------------rotration----------------------
        if uniform(0, 1) <= 0:
            disx = []
            disy = []
            disz = []
            for pot in cluster_dis.keys():
                disx.append(cluster_dis[pot][0])
                disy.append(cluster_dis[pot][1])
                disz.append(cluster_dis[pot][2])
            X = np.mean(disx)
            Y = np.mean(disy)
            Z = np.mean(disz)
            center = [X, Y, Z]

            dlist = []
            for i in range(1000):
                for pot_choice in cluster_dis.keys():
                    if pot_choice not in dlist and uniform(0, 1) < 0.333:
                        dlist.append(pot_choice)
                if len(dlist) >= 1:
                    break

            for pot0 in dlist:
                if 'Cu-' in pot0:
                    Theta = uniform(0, 2*pi)
                    intermediate = np.array(
                        cluster_dis[pot0]) - np.array(center)
                    Xarray = intermediate[0] * \
                        cos(Theta)-intermediate[1]*sin(Theta)
                    Yarray = intermediate[0] * \
                        sin(Theta)+intermediate[1]*cos(Theta)
                    Zarray = intermediate[2]
                    arr1_1 = [float(Xarray)+center[0], float(Yarray) +
                              center[1], float(cluster_dis[pot0][2])]
                    arr1 = np.array(arr1_1)
                    cluster_dis[pot0] = list(arr1_1)
                    allatom_dis[pot0] = list(arr1_1)

                    move_dis = [0, 0, 0]

                    for pot_cluster in cluster_gas.keys():
                        if pot0 == pot_cluster:
                            if not 'none' in cluster_gas[pot_cluster]:
                                for pot_gas in cluster_gas[pot_cluster].keys():
                                    intermediate_gas = np.array(
                                        cluster_gas[pot_cluster][pot_gas]) - np.array(center)
                                    X_gas = intermediate_gas[0] * \
                                        cos(Theta) - \
                                        intermediate_gas[1]*sin(Theta)
                                    Y_gas = intermediate_gas[0] * \
                                        sin(Theta) + \
                                        intermediate_gas[1]*cos(Theta)
                                    arr_gas = list([float(X_gas)+center[0]+move_dis[0], float(
                                        Y_gas)+center[1]+move_dis[1], float(cluster_gas[pot_cluster][pot_gas][2]+move_dis[2])])
                                    allatom_dis[pot_gas] = arr_gas
                                    cluster_gas[pot_cluster][pot_gas] = arr_gas
                    cluster_dis_new, cluster_gas_new = self.single_atom(
                        cluster_dis, cluster_gas)
                    allatom_dis.update(cluster_dis_new)
                    cluster_gas = cluster_gas_new
                    for pot in cluster_gas_new.keys():
                        if not 'none' in cluster_gas_new[pot]:
                            for pot_gas in cluster_gas_new[pot].keys():
                                allatom_dis[pot_gas] = cluster_gas_new[pot][pot_gas]

                else:
                    Theta = uniform(0, 2*pi)
                    intermediate = np.array(
                        cluster_dis[pot0]) - np.array(center)
                    Xarray = intermediate[0] * \
                        cos(Theta)-intermediate[1]*sin(Theta)
                    Yarray = intermediate[0] * \
                        sin(Theta)+intermediate[1]*cos(Theta)
                    Zarray = intermediate[2]
                    arr1_1 = [float(Xarray)+center[0], float(Yarray) +
                              center[1], cluster_dis[pot0][2]]
                    arr1 = np.array(arr1_1)
                    cluster_dis[pot0] = list(arr1_1)
                    allatom_dis[pot0] = list(arr1_1)

                    dis = 100
                    pot_1 = pot0
                    arr2_1 = np.array(cluster_dis[pot0])
                    move_dis = [0, 0, 0]
                    for pot_cluster in cluster_gas.keys():
                        if pot0 == pot_cluster:
                            if not 'none' in cluster_gas[pot_cluster]:
                                for pot_gas in cluster_gas[pot_cluster].keys():
                                    intermediate_gas = np.array(
                                        cluster_gas[pot_cluster][pot_gas]) - np.array(center)
                                    X_gas = intermediate_gas[0] * \
                                        cos(Theta) - \
                                        intermediate_gas[1]*sin(Theta)
                                    Y_gas = intermediate_gas[0] * \
                                        sin(Theta) + \
                                        intermediate_gas[1]*cos(Theta)
                                    arr_gas = list([float(X_gas)+center[0]+move_dis[0], float(
                                        Y_gas)+center[1]+move_dis[1], float(cluster_gas[pot_cluster][pot_gas][2]++move_dis[2])])
                                    allatom_dis[pot_gas] = arr_gas
                                    cluster_gas[pot_cluster][pot_gas] = arr_gas
                    cluster_dis_new, cluster_gas_new = self.single_atom(
                        cluster_dis, cluster_gas)
                    allatom_dis.update(cluster_dis_new)
                    cluster_gas = cluster_gas_new
                    for pot in cluster_gas_new.keys():
                        if not 'none' in cluster_gas_new[pot]:
                            for pot_gas in cluster_gas_new[pot].keys():
                                allatom_dis[pot_gas] = cluster_gas_new[pot][pot_gas]

        return allatom_dis, cluster_gas

# --------------add atom in cluster----------------------------------------------------
    def new_cluster_add(self, cluster, i, allatom1, cluster_gas):
        cluster_dis = copy.deepcopy(cluster)
        allatom = copy.deepcopy(allatom1)
        cluster_gas_dis = copy.deepcopy(cluster_gas)
        new_cluster = dict()
        disx = []
        disy = []
        disz = []
        for pot in cluster_dis.keys():
            disx.append(cluster_dis[pot][0])
            disy.append(cluster_dis[pot][1])
            disz.append(cluster_dis[pot][2])

        X = np.mean(disx)
        Y = np.mean(disy)
        Z = np.mean(disz)
        center = [X, Y, Z]
        for i in range(1000):
            if uniform(0, 1) <= 0.5:
                # --- we use r to search the position of O which is added----
                pot_all = []
                for pot in cluster_dis.keys():
                    pot_all.append(pot)
                pot_choice = choice(pot_all)
                dis_choice = cluster_dis[pot_choice]
                r = uniform(self.cutoff_distance3, self.cutoff_distance)
                arr_choice = list(np.array(dis_choice)-np.array(center))
                x_add = uniform(0, r)
                if arr_choice[0] > 0:
                    x_uniform = dis_choice[0] + x_add
                else:
                    x_uniform = dis_choice[0] - x_add

                d_y_z = sqrt(r**2-x_add**2)
                y_add = uniform(0, d_y_z)

                if arr_choice[1] > 0:
                    y_uniform = dis_choice[1] + y_add
                else:
                    y_uniform = dis_choice[1] - y_add

                d_z = sqrt(d_y_z**2-y_add**2)

                if uniform(0, 1) > 0.5:
                    z_uniform = dis_choice[0] + d_z
                else:
                    z_uniform = dis_choice[0] - d_z

                new_atom = [x_uniform, y_uniform, z_uniform]

            else:
                box = [min(disx)-2.5, max(disx)+2.5, min(disy)-2.5,
                       max(disy)+2.5, max(1.6, min(disz)), max(disz)+3]
                new_atom = [uniform(box[0], box[1]), uniform(
                    box[2], box[3]), uniform(box[4], box[5])]

            arr1 = np.array(new_atom)
            dis = 100
            for pot1 in cluster_dis.keys():
                arr2 = np.array(cluster_dis[pot1])
                dist = np.linalg.norm(arr1-arr2)
                if dist < dis:
                    dis = dist
                    arr2_1 = arr2
                    pot_1 = pot1

            if len(arr2_1) != 0:
                if dis > self.cutoff_distance or dis < self.cutoff_distance3:
                    if 'Cu-' in pot_1:
                        move_vector = arr1-arr2_1
                        move_vector = move_vector*self.cutoff_distance/dis
                        arr1 = arr2_1 + move_vector
                    else:
                        move_vector = arr1-arr2_1
                        move_vector = move_vector*self.cutoff_distance3/dis
                        arr1 = arr2_1 + move_vector

            new_cluster.update(cluster_dis)
            new_cluster['O_%s' % i] = list(arr1)
            allatom['O_%s' % i] = list(arr1)
            cluster_gas_dis['O_%s' % i] = 'none'
            # ----------------add cluster >= H_fix-------------
            # print(list(arr1)[2])
            if float(list(arr1)[2]) > max(2, min(disz)):
                break
        return allatom, cluster_gas_dis

# --------------add gas adoption new----------------------------------------------------------------------------------------------------------------
    def add_gas(self, support, cluster, t, allatom1, cluster_gas1):
        support_dis = support
        cluster_gas = copy.deepcopy(cluster_gas1)
        allatom = copy.deepcopy(allatom1)
        disx1 = []
        disy1 = []
        disz1 = []
        gas_new = dict()
        for pot in support_dis.keys():
            disx1.append(support_dis[pot][0])
            disy1.append(support_dis[pot][1])
            disz1.append(support_dis[pot][2])
        box = [min(disx1), max(disx1), min(disy1),
               max(disy1), min(disz1), max(disz1)]
        cluster_dis = cluster.copy()
        disx = []
        disy = []
        disz = []
        for pot in cluster_dis.keys():
            disx.append(cluster_dis[pot][0])
            disy.append(cluster_dis[pot][1])
            disz.append(cluster_dis[pot][2])
        box1 = [min(disx)-1.5, max(disx)+1.5, min(disy)-1.5,
                max(disy)+1.5, min(disz)+1, max(disz)+1]
        atom_gas_C = ['C', 'C']
        atom_gas_O = ['O2', 'O2', 'O2', 'O2', 'O2']
        x1 = [uniform(box[0], box1[0]), uniform(box[1], box1[1])]
        y1 = [uniform(box[2], box1[2]), uniform(box[3], box1[3])]
        z = float(uniform(box1[4], box1[5]))
        volume1 = (box[1]-box[0])*(box[3]-box[2])*(box1[5]-box1[4]) - \
            (box1[1]-box1[0])*(box1[3]-box1[2])*(box1[5]-box1[4])
        volume2 = (box[1]-box[0])*(box[3]-box[2])*(10-box1[5])
        pro_side = volume1/(volume1+volume2)
        for i in range(1000):
            if uniform(0, 1) <= pro_side:
                dis = [float(sample(x1, 1)[0]), float(sample(y1, 1)[0]), z]
            else:
                dis = [float(uniform(box[0], box[1])), float(
                    uniform(box[2], box[3])), float(uniform(box1[5], 15))]
            C_pro = uniform(0, 1)
            if C_pro <= self.ratio_CO:
                atom_gas = atom_gas_C
            else:
                atom_gas = atom_gas_O
            atom_i = sample(atom_gas, 1)
            arr1 = np.array(dis)
            ini = 100
            x = 0
            for pot in cluster_dis.keys():
                if 'Cu' in pot:
                    arr2_1 = np.array(cluster_dis[pot])
                    dist = np.linalg.norm(arr1-arr2_1)
                    if dist < ini:
                        ini = dist
                        arr2 = arr2_1
                        atom_cluster = pot
            for pot_get in cluster_gas.keys():
                if atom_cluster == pot_get:
                    if 'none' in cluster_gas[pot_get]:
                        x += 1
                        break

        # ----------------------if cluster_gas is all ad,continue ad--------
            x = 1
            if x == 1:
                if 'C' in atom_i:
                    move_vector = arr1-arr2
                    move_vector = move_vector / \
                        np.linalg.norm(move_vector)*self.cutoff_distance
                    arr1 = arr2+move_vector
                    O_C = arr1.copy()
                    move_vector = O_C-arr2
                    move_vector = move_vector / \
                        np.linalg.norm(move_vector)*self.cutoff_distance1
                    O_C = O_C+move_vector
                    arr1 = list(arr1)
                    O_C = list(O_C)
                    gas_new['C-%s' % t] = arr1
                    gas_new['O1-%s' % t] = O_C
                    allatom['C-%s' % t] = arr1
                    allatom['O1-%s' % t] = O_C
                    new_gas = {'C-%s' % t: arr1, 'O1-%s' % t: O_C}
                    cluster_gas[atom_cluster] = new_gas
                elif 'O1' in atom_i:
                    move_vector = arr1-arr2
                    move_vector = move_vector / \
                        np.linalg.norm(move_vector)*self.cutoff_distance
                    arr1 = arr2+move_vector
                    arr1 = list(arr1)
                    gas_new[atom_i[0]+'-%s' % t] = arr1
                    allatom[atom_i[0]+'-%s' % t] = arr1
                    cluster_gas[atom_cluster] = {atom_i[0]+'-%s' % t: arr1}
                elif 'O2' in atom_i:
                    move_vector = arr1-arr2
                    move_vector = move_vector / \
                        np.linalg.norm(move_vector)*self.cutoff_distance
                    arr1 = arr2+move_vector
                    O_O = arr1.copy()
                    move_vector = O_O-arr2
                    move_vector = move_vector / \
                        np.linalg.norm(move_vector)*self.cutoff_distance1
                    O_O = O_O+move_vector
                    arr1 = list(arr1)
                    O_O = list(O_O)
                    gas_new[atom_i[0]+'-1-%s' % t] = arr1
                    gas_new[atom_i[0]+'-2-%s' % t] = O_O
                    allatom[atom_i[0]+'-1-%s' % t] = arr1
                    allatom[atom_i[0]+'-2-%s' % t] = O_O
                    new_gas = {atom_i[0]+'-1-%s' %
                               t: arr1, atom_i[0]+'-2-%s' % t: O_O}
                    cluster_gas[atom_cluster] = new_gas
                break
        return allatom, cluster_gas

# --------------replace cluster adoption new----------------------------------------------------------------------------------------------------------------
    def re_cluster(self, cluster, allatom1, cluster_gas):
        cluster_dis = copy.deepcopy(cluster)
        allatom = copy.deepcopy(allatom1)
        cluster_gas_dis = copy.deepcopy(cluster_gas)
        new_cluster = dict()
        Cu_cluster = []
        O_cluster = []
        for pot in cluster_dis.keys():
            if 'Cu-' in pot:
                Cu_cluster.append(pot)
            elif 'O' in pot:
                O_cluster.append(pot)
        if len(O_cluster) >= 1:
            cluster_choice_O = sample(O_cluster, 1)[0]
            cluster_choice_Cu = sample(Cu_cluster, 1)[0]
            allatom[cluster_choice_Cu] = allatom1[cluster_choice_O]
            allatom[cluster_choice_O] = allatom1[cluster_choice_Cu]
            cluster_gas_dis[cluster_choice_Cu] = cluster_gas[cluster_choice_O]
            cluster_gas_dis[cluster_choice_O] = cluster_gas[cluster_choice_Cu]
        return allatom, cluster_gas_dis


# --------------replace gas adoption new----------------------------------------------------------------------------------------------------------------

    def re_gas(self, t, allatom1, cluster_gas1):
        cluster_gas = copy.deepcopy(cluster_gas1)
        allatom = copy.deepcopy(allatom1)
        cluster_gas_all = []
        x = 0
        for pot in cluster_gas.keys():
            if not 'none' in cluster_gas[pot]:
                cluster_gas_all.append(pot)
        if len(cluster_gas_all) >= 1:
            cluster_choice = sample(cluster_gas_all, 1)[0]
            arr1 = np.array(allatom[cluster_choice])
            dis_min, pot_1 = self.min_distance(
                arr1, cluster_choice, cluster_gas[cluster_choice])
            if 'O' in pot_1:
                number_gas = len(cluster_gas[cluster_choice])
                if number_gas >= 2:
                    arr1_choice = np.array(allatom[cluster_choice])
                    dis_min_min, pot_2 = self.min_distance(
                        arr1_choice, pot_1, cluster_gas[cluster_choice])
                    if 'O' in pot_2 and dis_min_min <= 1.5:
                        new_key = 'C-%s' % t
                        allatom.pop(pot_1)
                        allatom[new_key] = allatom1[pot_1]
                        cluster_gas[cluster_choice].pop(pot_1)
                        cluster_gas[cluster_choice][new_key] = allatom1[pot_1]
                    else:
                        x += 1
                else:
                    x += 1
            elif 'C-' in pot_1:
                new_key = 'O-%s' % t
                allatom.pop(pot_1)
                allatom[new_key] = allatom1[pot_1]
                cluster_gas[cluster_choice].pop(pot_1)
                cluster_gas[cluster_choice][new_key] = allatom1[pot_1]
            else:
                x += 1
        return allatom, cluster_gas, x

# --------------write new input.arc----------------------------------------------------------------------------------------------------------------------
    def new_input(self, atom_all, filename, newfilename):
        for line in range(1, 6):
            start = linecache.getline('./%s' % filename, line)
            with open('./%s' % newfilename, 'a') as fp1:
                fp1.write(start)
        a = 1
        for pot1 in atom_all.keys():
            atom = re.findall(r'[A-Za-z]', pot1)
            atom = ''.join(atom)
            with open('./%s' % newfilename, 'a') as fp1:
                fp1.write('%s  %15.9f%15.9f%15.9f CORE %s  %s %s 0.0000 %s\n' % (atom.ljust(3), float(atom_all[pot1][0]), float(
                    atom_all[pot1][1]), float(atom_all[pot1][2]), str('%4.0f' % (a)), atom.ljust(3), atom.ljust(3), str('%4.0f' % (a))))
            a += 1
        with open('./%s' % newfilename, 'a') as f:
            f.write('end\n')
            f.write('end')
        linecache.clearcache()

# -----------------------------same structure-------------------
    def is_same_structure(self, structure1, structure2):
        dis_all = []
        a = 0
        for pot1 in structure1.keys():
            for pot2 in structure2.keys():
                if pot1 == pot2:
                    arr1 = np.array(structure1[pot1])
                    arr2 = np.array(structure2[pot2])
                    distance = float(np.linalg.norm(arr1-arr2))
                    dis_all.append(distance)
        if sum(dis_all)/(len(structure1)+len(structure2)) < 0.015 and max(dis_all) < 0.4:
            a = 1
        else:
            a = 0
        return a

# --------------badstr-adsorption----------------------------------------------------------------------------------------------------------------------
    def bad1(self, gas, structure):
        gas_dis = copy.deepcopy(gas)
        structure_dis = copy.deepcopy(structure)
        a = 0
        for pot in gas_dis.keys():
            arr1 = np.array(gas_dis[pot])
            dis_min, pot_1 = self.min_distance(arr1, pot, structure_dis)
            if 'C-' in pot:
                if 'Cu-' in pot_1:
                    if dis_min >= 2.5 or dis_min < 1.6:
                        a += 1
                elif 'O' in pot_1:
                    if dis_min >= 2.1 or dis_min < 1.1:
                        a += 1
                elif 'C-' in pot_1:
                    if dis_min >= 2 or dis_min <= 1.6:
                        a += 1
            elif 'O' in pot:
                if 'Cu-' in pot_1:
                    if dis_min >= 2.5 or dis_min < 1.5:
                        a += 1
                elif 'C-' in pot_1:
                    if dis_min >= 2 or dis_min <= 1.1:
                        a += 1
                elif 'O' in pot_1:
                    if dis_min >= 2 or dis_min <= 1.1:
                        a += 1
        return a


# --------------badstr-cluster----------------------------------------------------------------------------------------------------------------------


    def bad2(self, structure):
        structure_dis = copy.deepcopy(structure)
        a = 0
        for pot in structure_dis.keys():
            arr1 = np.array(structure_dis[pot])
            dis_min, pot_1 = self.min_distance(arr1, pot, structure_dis)

            if 'Cu-' in pot:
                if 'O' in pot_1:
                    if dis_min >= 2.2 or dis_min < 1.3:
                        # print(dis_min)
                        a += 1

                elif 'Cu-' in pot_1:
                    if dis_min >= 2.8 or dis_min < 1.5:
                        # print(dis_min)
                        a += 1

                elif 'C-' in pot_1:
                    if dis_min >= 3 or dis_min < 1.6:
                        # print(dis_min)
                        a += 1

            elif 'O' in pot:
                if 'Cu-' in pot_1:
                    if dis_min >= 2.5 or dis_min < 1.3:
                        # print(dis_min)
                        a += 1

                elif 'O' in pot_1:
                    if dis_min >= 3 or dis_min < 0.9:
                        # print(dis_min)
                        a += 1

                elif 'C-' in pot_1:
                    if dis_min >= 3 or dis_min < 1:
                        # print(dis_min)
                        a += 1

            elif 'C-' in pot:
                if 'Cu-' in pot_1:
                    if dis_min >= 3 or dis_min < 1.6:
                        a += 1

                elif 'O' in pot_1:
                    if dis_min >= 3 or dis_min < 1:
                        a += 1

        return a

# --------------badstr-gas and gas(bond length)----------------------------------------------------------------------------------
    def bad3(self, gas):
        gas_dis = copy.deepcopy(gas)
        if len(gas_dis) > 1:
            dis = 100
            a = 0
            for pot in gas_dis.keys():
                arr1 = np.array(gas_dis[pot])
                dis_min, pot_1 = self.min_distance(arr1, pot, gas_dis)
                if 'O' in pot:
                    if 'C-' in pot_1:
                        if dis_min < 1.1:
                            a += 1
                elif 'O' in pot:
                    if 'O' in pot_1:
                        if dis_min < 1.1:
                            a += 1
                elif 'C-' in pot:
                    if 'C-' in pot_1:
                        if dis_min < 2:
                            a += 1

                elif 'C-' in pot:
                    if 'O' in pot_1:
                        if dis_min < 1.1:
                            a += 1
        else:
            a = 0
        return a

# --------------badstr-gas and gas(distance)----------------------------------------------------------------------------------
    def bad4(self, gas):
        gas_dis = gas.copy()
        a = 0
        if len(gas_dis) >= 2:
            for pot in gas_dis.keys():
                dis1 = 100
                arr1 = np.array(gas_dis[pot])
                for pot1 in gas_dis.keys():
                    if pot == pot1:
                        continue
                    else:
                        arr2 = np.array(gas_dis[pot1])
                        dist = np.linalg.norm(arr1-arr2)
                        if dist < dis1:
                            dis1 = dist
                            pot_0 = pot
                            pot_1 = pot1

                dis2 = 100
                for pot11 in gas_dis.keys():
                    if pot_0 == pot11:
                        continue
                    elif pot11 == pot_1:
                        continue
                    else:
                        arr2 = np.array(gas_dis[pot11])
                        dist = np.linalg.norm(arr1-arr2)
                        if dist < dis2:
                            dis2 = dist
                            pot_2 = pot11
                if 'O' in pot_0:
                    if 'C-' in pot_2:
                        if dis2 <= 2.3:
                            a += 1

                elif 'O' in pot_0:
                    if 'O' in pot_2:
                        if dis2 <= 2.3:
                            a += 1

                elif 'C-' in pot_0:
                    if 'C-' in pot_2:
                        if dis2 <= 2.3:
                            a += 1

                elif 'O' in pot_0:
                    if 'O' in pot_2:
                        if dis2 <= 2.3:
                            a += 1

        return a

# --------------badstr-single gas or CO3(more) not gas----------------------------------------------------------------------------------
    def bad5(self, gas):
        gas_dis = gas.copy()
        a = 0
        if len(gas_dis) > 2:
            for pot in gas_dis.keys():
                dis = 17
                arr1 = np.array(gas_dis[pot])
                for pot1 in gas_dis.keys():
                    if pot != pot1:
                        arr2 = np.array(gas_dis[pot1])
                        dist = np.linalg.norm(arr1-arr2)
                        if dist < dis:
                            dis = dist
                            pot_1 = pot1
                        else:
                            continue
                dis = 17
                for pot1 in gas_dis.keys():
                    if pot != pot1 and pot1 != pot_1:
                        arr2 = np.array(gas_dis[pot1])
                        dist = np.linalg.norm(arr1-arr2)
                        if dist < dis:
                            dis = dist
                            pot_2 = pot1
                dis = 1.7
                for pot1 in gas_dis.keys():
                    if pot != pot1 and pot1 != pot_1 and pot1 != pot_2:
                        arr2 = np.array(gas_dis[pot1])
                        dist = np.linalg.norm(arr1-arr2)
                        if dist < dis:
                            dis = dist
                            pot_2 = pot1

                if dis <= 1.6:
                    a += 1

        return a

# --------------del-gas----------------------------------------------------------------------------------------------------------------------
    def del_gas(self, gas, allatom, cluster_gas):
        gas_dis = copy.deepcopy(gas)
        allatom_dis = copy.deepcopy(allatom)
        cluster_gas_dis = copy.deepcopy(cluster_gas)
        choose = []
        for pot in cluster_gas_dis.keys():
            if not 'none' in cluster_gas_dis[pot]:
                choose.append(pot)
        if len(choose) >= 1:
            atom_cluster = sample(choose, 1)[0]
            for pot_gas in cluster_gas_dis[atom_cluster].keys():
                del allatom_dis[pot_gas]
            if len(cluster_gas_dis[atom_cluster]) == 0:
                cluster_gas_dis[atom_cluster] = 'none'
        return allatom_dis, cluster_gas_dis

# --------------del-cluster---------------------------------------------------------------------------------------------------
    def del_cluster(self, cluster, allatom, cluster_gas):
        cluster_dis = copy.deepcopy(cluster)
        allatom_dis = copy.deepcopy(allatom)
        cluster_gas_dis = copy.deepcopy(cluster_gas)
        potO = []
        potCu = []
        for pot in cluster_dis.keys():
            if 'O-' in pot:
                potO.append(pot)
            elif 'Cu' in pot:
                potCu.append(pot)
        atom = sample(potO, 1)[0]
        atom_dis = [cluster[atom][0], cluster[atom][1], cluster[atom][2]]
        del cluster_dis[atom]
        del allatom_dis[atom]
        if not 'none' in cluster_gas_dis[atom]:
            for pot_gas in cluster_gas_dis[atom]:
                del allatom_dis[pot_gas]
        del cluster_gas_dis[atom]

        cluster_gas_dis[atom] = 'none'

        return allatom_dis, cluster_gas_dis

# -------------- if not CO2 had get---------------------------------------------------------------------------------------------------
    def CO2_bave(self, allatom):
        allatom_dis1 = copy.deepcopy(allatom)
        atom_CO2 = []
        x = 0
        for pot in allatom.keys():
            allatom_dis = copy.deepcopy(allatom_dis1)
            if 'C_' in pot or 'C-' in pot:
                arr1 = np.array(allatom_dis[pot])
                dis_min1, pot_1 = self.min_distance(arr1, pot, allatom_dis)
                if 'O-' in pot_1 and dis_min1 <= 1.7:
                    del allatom_dis[pot_1]
                    dis_min2, pot_2 = self.min_distance(arr1, pot, allatom_dis)
                    if 'O-' in pot_2 and dis_min2 <= 1.7:
                        del allatom_dis[pot_2]
                        dis_min3, pot_3 = self.min_distance(
                            arr1, pot, allatom_dis)

                        if dis_min3 >= 2 or 'Cu' in pot_3:
                            x += 1
                            atom_CO2.append(pot)
                            atom_CO2.append(pot_1)
                            atom_CO2.append(pot_2)
                        elif 'O' in pot_3 and dis_min3 <= 1.7:
                            x += 1
                            atom_CO2.append(pot)
                            atom_CO2.append(pot_1)
                            atom_CO2.append(pot_2)
                    elif dis_min2 > 2.5:
                        x += 1
                        atom_CO2.append(pot)
                        atom_CO2.append(pot_1)

        if x != 0:
            for pot_del in atom_CO2:
                del allatom_dis1[pot_del]

        return x, allatom_dis1

# ------------- del sur gas ---------------------------------------------------------------------------------------------------
    def del_sur_gas(self, cluster, gas, allatom):
        gas_dis = copy.deepcopy(gas)
        gas_dis_copy = copy.deepcopy(gas)
        cluster_dis = copy.deepcopy(cluster)
        structure_dis = copy.deepcopy(allatom)

        if len(gas_dis) >= 2:
            for pot in gas_dis.keys():
                dis = 100
                arr1 = np.array(gas_dis[pot])
                for pot1 in structure_dis.keys():
                    if pot == pot1:
                        continue
                    else:
                        arr2 = np.array(structure_dis[pot1])
                        dist = np.linalg.norm(arr1-arr2)
                        if dist < dis:
                            dis = dist
                            pot_0 = pot
                            pot_1 = pot1
                        else:
                            continue
                if 'C-' in pot_0:
                    if 'Cu' in pot_1:
                        if dis >= 2.7:
                            del gas_dis_copy[pot_0]
                    elif 'O' in pot_1:
                        if dis >= 2.5:
                            del gas_dis_copy[pot_0]

                elif 'O' in pot_0:
                    if dis >= 5:
                        del gas_dis_copy[pot_0]
                    elif dis <= 5 and dis >= 2.5:
                        if len(gas_dis) > 1:
                            dis = 100
                            arr1 = np.array(gas_dis[pot_0])
                            for pot1 in gas_dis.keys():
                                if pot_0 == pot1:
                                    continue
                                else:
                                    arr2 = np.array(gas_dis[pot1])
                                    dist = np.linalg.norm(arr1-arr2)
                                    if dist < dis:
                                        dis = dist
                                        pot_1_1 = pot1
                            if dis <= 1.6 and dis > 1.2:
                                dis_ini = 100
                                arr1 = np.array(gas_dis[pot_1_1])
                                for pot1 in cluster_dis.keys():
                                    arr2 = np.array(cluster_dis[pot1])
                                    dist = np.linalg.norm(arr1-arr2)
                                    if dist < dis_ini:
                                        dis_ini = dist
                                        pot_1 = pot1
                                if dis_ini >= 2.7:
                                    del gas_dis_copy[pot_0]

            return structure_dis


# -------------------------badall----------------------------

    def badall(self, atom_all, H):
        support_atom_post, Height_post = self.support_atom(atom_all, H)
        cluster_atom_post = self.cluster_atom(atom_all, Height_post)
        gas_move = self.gas_atom(
            atom_all, support_atom_post, cluster_atom_post)
        cluster_and_gas = dict()
        cluster_and_gas.update(cluster_atom_post)
        cluster_and_gas.update(gas_move)
        if gas_move is not None:
            if len(gas_move) > 2:
                a1 = self.bad1(gas_move, atom_all)
                a2 = self.bad2(cluster_and_gas)
                a3 = self.bad3(gas_move)
                a4 = self.bad4(gas_move)
            else:
                a1 = 0
                a2 = self.bad2(cluster_and_gas)
                a3 = 0
                a4 = 0
        else:
            a1 = 0
            a2 = self.bad2(cluster_and_gas)
            a3 = 0
            a4 = 0
        print(a1, a2, a3, a4)
        a = a1+a2+a3+a4
        return a

# -------------------G with P-------------------------------------
    def G_with_P(self, energy, ratio):
        G_P = energy + 8.314*self.temperature * \
            log((self.pressure*ratio/self.standard_pressure))/(1000*96.4916)
        return G_P

# ------------- energy correction ---------------------------------------------------------------------------------------------------
    def energy_correction(self, num_C, num_O, num_CO2, num_cluster, energy):
        if self.ratio_CO != 0 and self.ratio_CO != 1:
            CO_energy = self.G_with_P(
                self.potential_CO, self.ratio_CO) + self.adjust_value
            O2_energy = self.G_with_P(
                self.potential_O2, (1-self.ratio_CO)) + self.adjust_value
        elif self.ratio_CO == 0:
            CO_energy = self.potential_CO + self.adjust_value
            O2_energy = self.G_with_P(
                self.potential_O2, (1-self.ratio_CO)) + self.adjust_value
        elif self.ratio_CO == 1:
            O2_energy = self.potential_O2 + self.adjust_value
            CO_energy = self.G_with_P(
                self.potential_CO, self.ratio_CO) + self.adjust_value
        Energy = energy + (self.all_CO-num_C-num_CO2)*CO_energy+(self.all_O2-num_O -
                                                                 num_CO2)*O2_energy - (num_cluster-3)*O2_energy + num_CO2*self.potential_CO2
        return Energy

# ------------- MC energy to get true str ---------------------------------------------------------------------------------------------------
    def energy_estimate(self, energy, ini_energy):
        a = 1
        if energy <= ini_energy:
            a = 0
        else:
            probability = exp(-1/(self.ConstantR*self.temperature)
                              * (energy - ini_energy))
            if probability > uniform(0, 1):
                a = 0
            else:
                a = 1
        return a

# --------------------error data about cluster-----------
    def err_data(self, cluster):
        x_data = []
        y_data = []
        z_data = []
        for pot in cluster.keys():
            x = cluster[pot][0]
            y = cluster[pot][1]
            z = cluster[pot][2]
            x_data.append(x)
            y_data.append(y)
            z_data.append(z)
        x_mean = np.mean(x_data)
        x_std = np.std(x_data, ddof=1)
        y_mean = np.mean(y_data)
        y_std = np.std(y_data, ddof=1)
        z_mean = np.mean(z_data)
        z_std = np.std(z_data, ddof=1)
        a = 0
        for pot in cluster.keys():
            x = cluster[pot][0]
            y = cluster[pot][1]
            z = cluster[pot][2]
            if x < (x_mean-2*x_std) or x > (x_mean+2*x_std) or y < (y_mean-2*y_std) or y > (y_mean+2*y_std) or z < (z_mean-2*z_std) or z > (z_mean+2*z_std):
                a += 1
                break
            else:
                a += 0
        return a

# --------------cluster-dis-rational----------------------------------------------
    def move_all_cluster(self, cluster, support):
        cluster_dis = cluster.copy()
        support_dis = support.copy()
        disx = []
        disy = []
        disz = []
        for pot in cluster_dis.keys():
            disx.append(cluster_dis[pot][0])
            disy.append(cluster_dis[pot][1])
            disz.append(cluster_dis[pot][2])
        X = np.mean(disx)
        Y = np.mean(disy)
        Z = np.mean(disz)
        center = [X, Y, 0]
        Z_min = min(disz)
        disx_sup = []
        disy_sup = []
        disz_sup = []
        for pot1 in support_dis.keys():
            disx_sup.append(support_dis[pot1][0])
            disy_sup.append(support_dis[pot1][1])
            disz_sup.append(support_dis[pot1][2])
        X_sup = np.mean(disx_sup)
        Y_sup = np.mean(disy_sup)
        Z_sup = np.mean(disz_sup)
        Z_max = max(disz_sup)
        # -------------z no change---------------
        center_sup = [X_sup, Y_sup, 0]
        arr_clu_to_sup = np.array(center_sup) - np.array(center)

        if abs(Z_min-Z_max) >= 2:
            arr_clu_to_sup[2] += 1.6-Z_min+Z_max
            for pot in cluster_dis.keys():
                cluster_dis[pot] = list(
                    np.array(cluster_dis[pot]) + arr_clu_to_sup)

        return cluster_dis

# --------------cluster-if-or-not-rational-----------------------------------------------
    def cluster_ration(self, cluster, Height):
        # cluster_dis = cluster.copy()
        cluster_dis = copy.deepcopy(cluster)

        x_final = 0
        disx = []
        disy = []
        disz = []
        dis_all = []
        for pot in cluster_dis.keys():
            dis_all.append([cluster[pot][0], cluster[pot][1], cluster[pot][2]])
        for i in range(1):
            for pot in cluster_dis.keys():
               #     if 'Cu' in pot and uniform(0,1) <= 0.9:
                if 'Cu-' in pot:
                    # if uniform(0,1) <= 1:
                    disx.append(cluster_dis[pot][0])
                    disy.append(cluster_dis[pot][1])
                    disz.append(cluster_dis[pot][2])
            center = [np.mean(disx), np.mean(disy), np.mean(disz)]
            x_std = np.std(disx, ddof=1)
            y_std = np.std(disy, ddof=1)
            z_std = np.std(disz, ddof=1)
            for i in range(len(cluster_dis)):
                x = dis_all[i][0]
                y = dis_all[i][1]
                z = dis_all[i][2]
                if x < (center[0]-3*x_std) or x > (center[0]+3*x_std) or y < (center[1]-3*y_std) or y > (center[1]+3*y_std) or z < (center[2]-3*z_std) or z > (center[2]+3*z_std) or z < Height or z > Height+6:
                    # if x <= center[0]-4.5 or x >= center[0]+5 or y <= center[0]-5 or y >= center[0]+5 or z <= center[0]-5 or x >= center[0]+5 or z<Height or z>Height+3:
                    x_final += 1
                    break
            if x_final == 1:
                break
        return x_final

# -------------------------get CO/O2 from gas (list)------------------------
    def get_CO_O2(self, gas):
        gas_dis = gas
        gas_copy = gas.copy()
        gas_copy = copy.deepcopy(gas)
        gas_all_list = []
        if len(gas_dis) >= 2:
            for pot1 in gas_dis.keys():
                arr1 = np.array(gas_dis[pot1])
                dis_min, pot_2 = self.min_distance(arr1, pot1, gas)
                if dis_min <= self.cutoff_distance:
                    for pot in gas_copy.keys():
                        if pot1 == pot:
                            gas_two = dict()
                            gas_two[pot1] = gas_dis[pot1]
                            del gas_copy[pot1]
                            a = 1
                            for pot_post in gas_copy.keys():
                                if pot_post == pot_2:
                                    break
                                else:
                                    a += 1
                            if a <= len(gas_copy):
                                gas_two[pot_2] = gas_dis[pot_2]
                                del gas_copy[pot_2]
                            gas_all_list.append(gas_two)
                            break

            if len(gas_copy) > 0:
                for pot1 in gas_copy.keys():
                    gas_single = dict()
                    gas_single[pot1] = gas_copy[pot1]
                    gas_all_list.append(gas_single)

        elif len(gas_dis) == 1:
            for pot1 in gas_copy.keys():
                gas_single = dict()
                gas_single[pot1] = gas_copy[pot1]
                gas_all_list.append(gas_single)
        else:
            gas_all_list.append('none')
        return gas_all_list

# --------------------------min distance(same stucture)-----------------------------------------------
    def min_distance(self, arr, pot1, struture):
        ini = 100
        arr1 = arr
        for pot2 in struture.keys():
            if not pot1 == pot2:
                arr2 = np.array(struture[pot2])
                dist = np.linalg.norm(arr1-arr2)
                if dist < ini:
                    ini = dist
                    pot_2 = pot2
        return ini, pot_2

# --------------------------cluster which adsorp gas (cluster:gas/none(dict))-----------------------------------------------
    def cluster_adsorp(self, structure, H):
        structure = copy.deepcopy(structure)
        support, Height = self.support_atom(structure, H)
        cluster = self.cluster_atom(structure, Height)
        gas = self.gas_atom(structure, support, cluster)
        gas_all_list = self.get_CO_O2(gas)
        cluster_ad = list()
        cluster_and_gas = dict()
        gas_no_choose = []
        number_cluster_ad = 0
        if not 'none' in gas_all_list:
            for i, pot_all in enumerate(gas_all_list):
                number_cluster_ad = 0
                for pot1 in pot_all.keys():
                    arr1 = np.array(pot_all[pot1])
                    dis_min, pot_2 = self.min_distance(arr1, pot1, cluster)
                    if dis_min <= 2.3 and dis_min >= 1.5:
                        cluster_ad.append(pot_2)
                        cluster_and_gas[pot_2] = pot_all
                        break
                    else:
                        number_cluster_ad += 1
                if number_cluster_ad == len(pot_all):
                    gas_no_choose.append(pot_all)
            for pot2 in cluster.keys():
                x = 1
                for pot_ad in cluster_and_gas.keys():
                    if not pot2 == pot_ad:
                        x += 1
                    else:
                        break
                if x > len(cluster_and_gas):
                    cluster_and_gas[pot2] = 'none'
        else:
            for pot2 in cluster.keys():
                cluster_and_gas[pot2] = 'none'
        return cluster_and_gas, gas_no_choose

# -----------------------------gas2cluster-------------------
    def gas2cluster(self, structure, cluster_and_gas):
        cluster_gas = copy.deepcopy(cluster_and_gas)
        structure_dis = copy.deepcopy(structure)
        for pot in cluster_gas.keys():
            if not 'none' in cluster_gas[pot]:
                arr_cluster = np.array(structure_dis[pot])
                dis_min, pot_2 = self.min_distance(
                    arr_cluster, pot, cluster_gas[pot])
                if dis_min >= 2 or dis_min <= 1.5:
                    move_vector = arr_cluster-np.array(cluster_gas[pot][pot_2])
                    move_vector_post = move_vector/dis_min * \
                        self.cutoff_distance+np.array(cluster_gas[pot][pot_2])
                    move_distance = move_vector_post - move_vector
                    for pot_gas in cluster_gas[pot].keys():
                        arr_gas = np.array(cluster_gas[pot][pot_gas])
                        arr1 = arr_gas+move_distance
                        cluster_gas[pot][pot_gas] = list(arr1)
                        structure[pot_gas] = list(arr1)
        return structure_dis, cluster_gas

# ------------------------------single atom-----------------------
    def single_atom(self, cluster, cluster_gas):
        cluster_dis = copy.deepcopy(cluster)
        cluster_gas_dis = copy.deepcopy(cluster_gas)
        for pot in cluster_dis.keys():
            arr_cluster = np.array(cluster_dis[pot])
            dis_min, pot_2 = self.min_distance(arr_cluster, pot, cluster)
            if dis_min >= self.cutoff_distance2:
                move_vector = arr_cluster-np.array(cluster_dis[pot_2])
                move_vector_post = move_vector / \
                    dis_min*(self.cutoff_distance2-0.2)
                cluster_dis[pot] = list(
                    move_vector_post + np.array(cluster_dis[pot_2]))
                move_distance = move_vector_post - move_vector
                if not 'none' in cluster_gas_dis[pot]:
                    for pot_gas in cluster_gas_dis[pot].keys():
                        arr_gas = np.array(cluster_gas_dis[pot][pot_gas])
                        arr1 = arr_gas+move_distance
                        cluster_gas_dis[pot][pot_gas] = list(arr1)

        return cluster_dis, cluster_gas_dis

    def CO_in_cluster(self, cluster, allatom):
        cluster_dis = copy.deepcopy(cluster)
        allatom_dis = copy.deepcopy(allatom)
        for pot in cluster_dis.keys():
            if 'C-' in pot:
                del allatom_dis[pot]
        return allatom_dis

    # ---------------atom out of the lattice, we need to modify atom----------
    def isolated_atom(self, atom, length, width):
        atom_dis = copy.deepcopy(atom)
        for pot in atom.keys():
            pot_dis = np.array(atom[pot])
            dis_min, pot_min = self.min_distance(pot_dis, pot, atom)
            if dis_min >= 7:
                x = float(list(pot_dis-np.array(atom[pot_min]))[0])
                y = float(list(pot_dis-np.array(atom[pot_min]))[1])
                z = float(list(pot_dis-np.array(atom[pot_min]))[2])
                if x >= 6:
                    pot_dis[0] = pot_dis[0] - length
                elif x <= -6:
                    pot_dis[0] = pot_dis[0] + length
                elif y >= 6:
                    pot_dis[1] = pot_dis[1] - width
                elif y <= -6:
                    pot_dis[1] = pot_dis[1] + width
                new_atom_pot_dis = np.array(
                    [pot_dis[0], pot_dis[1], pot_dis[2]])
                dis_min, pot_min = self.min_distance(
                    new_atom_pot_dis, pot, atom)
                print(dis_min, pot_min, new_atom_pot_dis)
                if dis_min <= 3:
                    atom_dis[pot] = list(new_atom_pot_dis)
        return atom_dis


def random_2part(number0, number1, number2, number3):
    if uniform(0, 1) <= (number1-number0)/(number1-number0+number3-number2):
        x = uniform(number0, number1)
    else:
        x = uniform(number2, number3)
    return x


def run(
        calc,
        file_meta='meta.traj',
        potential_O2=-1.071316E+01/2,
        potential_CO=-15.78612,
        potential_CO2=-24.24603,
        temperature=300,
        pressure=1e5,
        cycle=1000,
        stop_criterion=300,
        H_fix=1.6,
        ratio_CO=0.5):
    path = os.getcwd()
    # os.system('rm allarc; rm GCMCtraj;rm step-all;rm alllammps;rm energy-all;rm lammps-choice;rm alltrain;rm -r poscar-all;rm -r all-structure')
    # os.system('cp input.arc beststr.arc;cp input.arc allstruct.arc;cp input.arc allstr.arc; cp input.arc allstruct_gas.arc')
    os.system(
        'rm allposcar.traj;rm POSCAR-all;rm GCMCtraj;rm step-all;rm bestPOSCAR;rm best-POSCAR;rm time_all')
    path1 = path + '/best-1.arc'
    # pop = Population()
    pop = Population(ratio_CO=ratio_CO,
                      potential_O2=potential_O2,
                      potential_CO=potential_CO,
                      potential_CO2=potential_CO2,
                      temperature=temperature,
                      pressure=pressure)
    Path_arc = path + '/str'
    Path_arc2 = path + '/str2'
    Path_arc3 = path + '/str3'
    Path_lammps = path + '/lammps-str'
    if os.path.exists('./opt'):
        os.system('rm -r opt')
    os.mkdir('./opt')

# --------------ini --------------
    x = 0
    # allatom = pop.atom(path,'input.arc')
    all_times = 0
    # cycle = 20000
    # stop_criterion = 10000
    # H_fix = 1.6  # fix height
    # --------------the type of atom--------------
    # atomtype = ['Cu', 'Ce', 'O', 'C']
    # # -----------------the device is cpu or gpu( cpu is default)---------------
    # device = 'cpu'
    # # -------------------------pbc([1,1,1] is default)---------------------
    # period = [1, 1, 1]
    # # ---------------nn file('EANN_PES_DOUBLE.pt' is default)----------------------------
    # nn = 'EANN_PES_DOUBLE.pt'
    traj = Trajectory(file_meta)
    atoms = traj[0]
    #cell1 = ase.io.vasp.read_vasp("POSCAR")
    #atoms = Atoms(cell1)
    atoms.calc = calc
    dyn = LBFGS(atoms, trajectory='atom2.traj')
    dyn.run(fmax=0.1, steps=100)
    traj = Trajectory('atom2.traj')
    atoms = traj[-1]
    ase.io.write('best-POSCAR', atoms, format='vasp', vasp5='True')
    allatom_gas = POSCAR_To_atom1(path, 'best-POSCAR')
    energy_ini = atoms.get_potential_energy()
    if os.path.exists('./allstr.traj'):
        os.remove('./allstr.traj')
    traj = Trajectory('allstr.traj', 'a')
    traj.write(atoms)
    #stop_time = 0
    # energy_ini = -450.41673138151174   #ini_energy
    #etol = 0.1
    #ftol = 0.4
    #opt_max = 200
    # bfgs_max = 2         #-opt_max_time
    os.system('python dire2cart.py POSCAR;mv POSCAR_C POSCAR')
    allatom = POSCAR_To_atom1(path, 'POSCAR')
    support_atom, Height = pop.support_atom(allatom, H_fix)
    cluster_atom = pop.cluster_atom(allatom, Height)
    gas_atom_ini = pop.gas_atom(allatom, support_atom, cluster_atom)
    cluster_gas, gas_no_choose = pop.cluster_adsorp(allatom, H_fix)
    if len(gas_no_choose) > 0:
        allatom_copy = copy.deepcopy(allatom)
        for pot in gas_no_choose:
            for pot_gas in pot:
                del allatom_copy[pot_gas]
        allatom = allatom_copy

    num_atom_c = 0
    if len(gas_atom_ini) != 0:
        all_gas_num = len(gas_atom_ini)
        for pot in gas_atom_ini.keys():
            if 'C' in pot:
                num_atom_c += 1
        num_atom_O = all_gas_num - 2*num_atom_c
    else:
        num_atom_O = 0
    ini_energy = pop.energy_correction(
        num_atom_c, num_atom_O, 0, len(cluster_atom), energy_ini)
    stop_time = 0
    #best_energy = pop.energy_correction(num_atom_c, num_atom_O, 0, len(cluster_atom), energy_ini)
    #bestPOSCAR = -810.5575263
    #bestPOSCAR = pop.energy_correction(num_atom_c, num_atom_O, 0, len(cluster_atom), energy_ini)
    bestPOSCAR = ini_energy
# --------------1000 times MC--------------
    for cycle_time in range(cycle):
        x = 0
        if all_times <= cycle:
            times = all_times
            all_times += 1
        else:
            print('all have finished')
            break
# ------------------------ ini energy  -------------------------------
        traj = Trajectory('meta.traj')
        atoms = traj[randint(0, len(traj)-1)]
        ase.io.write('POSCAR', atoms, format='vasp', vasp5='True')
        os.system('python dire2cart.py POSCAR;mv POSCAR_C POSCAR')
        allatom = POSCAR_To_atom1(path, 'POSCAR')
        # allatom = pop.atom(path,'input.arc')
        support_atom, Height = pop.support_atom(allatom, H_fix)
        cluster_atom = pop.cluster_atom(allatom, Height)
        gas_atom_ini = pop.gas_atom(allatom, support_atom, cluster_atom)
        cluster_gas, gas_no_choose = pop.cluster_adsorp(allatom, H_fix)
# ----------------------------preprocessing-------------------------------
        length = 15.3045997620
        width = 13.2541721886
        atom_no_sup = copy.deepcopy(cluster_atom)
        atom_no_sup.update(gas_atom_ini)
        allatom_nosup_new = pop.isolated_atom(atom_no_sup, length, width)
        allatom_new1 = support_atom
        allatom_new1.update(allatom_nosup_new)
        allatom = allatom_new1

        support_atom, Height = pop.support_atom(allatom, H_fix)
        cluster_atom = pop.cluster_atom(allatom, Height)
        gas_atom_ini = pop.gas_atom(allatom, support_atom, cluster_atom)
        cluster_gas, gas_no_choose = pop.cluster_adsorp(allatom, H_fix)
        if len(gas_no_choose) > 0:
            allatom_copy = copy.deepcopy(allatom)
            for pot in gas_no_choose:
                for pot_gas in pot:
                    del allatom_copy[pot_gas]
            allatom = allatom_copy

# ----------------------------preprocessing-------------------------------
        atom_preprocessing, cluster_gas_preprocessing = pop.gas2cluster(
            allatom, cluster_gas)
        allatom = atom_preprocessing
        cluster_gas = cluster_gas_preprocessing

        support_atom, Height = pop.support_atom(allatom, H_fix)
        cluster_atom = pop.cluster_atom(allatom, Height)
        gas_atom_ini = pop.gas_atom(allatom, support_atom, cluster_atom)

        # ------------sur-atom----------------------
        atom_all_new = pop.del_sur_gas(cluster_atom, gas_atom_ini, allatom)
        atom_all = atom_all_new
        support_atom, Height = pop.support_atom(allatom, H_fix)
        cluster_atom = pop.cluster_atom(allatom, Height)
        gas_atom_ini = pop.gas_atom(allatom, support_atom, cluster_atom)
        # ---------------------------------------------
        cluster_gas, gas_no_choose = pop.cluster_adsorp(allatom, H_fix)
        if len(gas_no_choose) > 0:
            allatom_copy = copy.deepcopy(allatom)
            for pot in gas_no_choose:
                for pot_gas in pot:
                    del allatom_copy[pot_gas]
                    allatom = allatom_copy

        # atom_all_new = pop.del_sur_gas(cluster_atom,gas_atom_ini,allatom)

        num_atom_c = 0
        if len(gas_atom_ini) != 0:
            all_gas_num = len(gas_atom_ini)
            for pot in gas_atom_ini.keys():
                if 'C' in pot:
                    num_atom_c += 1
            num_atom_O = all_gas_num - 2*num_atom_c
        else:
            num_atom_O = 0

        with open('GCMCtraj', 'a+') as fp:
            fp.write('%s:\t%s\n' % (times, ini_energy))
        num_change = 0
# --------------------- input -> input.arc------------------
        os.system('cp POSCAR POSCAR-old')
# ------------------start MC---------------
        cluster_add = cluster_atom
        number1 = uniform(0, 1)
        #number1 = 1
        #probality = 1/3*(1-(len(cluster_add)-8)/9)
        probality = 1/3*0.5
        if len(cluster_add) < 4:
            if uniform(0, 1) < 0.25:
                number1 = 0 or 1
            else:
                number1 = uniform(1/3+0.01, 1)
# ###------------------get cluster add-----------------
        if number1 < probality:
            for i in range(2000):
                t = 1000 + times * 5
                allatom_copy = copy.deepcopy(allatom)
                cluster_gas_copy = copy.deepcopy(cluster_gas)
                atom_all, cluster_gas_new = pop.new_cluster_add(
                    cluster_atom, t, allatom, cluster_gas)
                a_bad = pop.badall(atom_all, H_fix)
                if a_bad == 0:
                    allatom = atom_all
                    cluster_gas = cluster_gas_new
                    with open('step-all', 'a') as fp:
                        fp.write(
                            '------%s-add-cluster-%s is finishing-------\n' % (times, i))
                    print('------%s-add-cluster-%s is finishing-------' %
                          (times, i))
                    break
                else:
                    if i < 1999:
                        allatom = allatom_copy
                        cluster_gas = cluster_gas_copy
                    else:
                        allatom = allatom_copy
                        cluster_gas = cluster_gas_copy
                        print('add-cluster-bad')
                        x += 1
                        break

# ###------------------get cluster del-----------------
        elif number1 >= probality and number1 <= 1/3:
            for i in range(2000):
                allatom_copy = copy.deepcopy(allatom)
                cluster_gas_copy = copy.deepcopy(cluster_gas)
                atom_all, cluster_gas_post = pop.del_cluster(
                    cluster_atom, allatom, cluster_gas)
                a_bad = pop.badall(atom_all, H_fix)
                if a_bad == 0:
                    x += 0
                    allatom = atom_all
                    cluster_gas = cluster_gas_post
                    with open('step-all', 'a') as fp:
                        fp.write(
                            '------%s-del cluster-%s is finishing-------\n' % (times, i))
                    print('------%s-del cluster-%s is finishing-------' %
                          (times, i))
                    break
                else:
                    if i < 1999:
                        allatom = allatom_copy
                        cluster_gas = cluster_gas_copy
                    else:
                        allatom = allatom_copy
                        cluster_gas = cluster_gas_copy
                        print('del-cluster-bad')
                        x += 1
                        break

# ###------------------get cluster move and rotation-----------------
        elif number1 > 1/3 and number1 <= 2/3:
            for i in range(2000):
                allatom_copy = copy.deepcopy(allatom)
                cluster_gas_copy = copy.deepcopy(cluster_gas)
                atom_all, cluster_gas_new = pop.new_cluster_move(
                    cluster_atom, Height, allatom, cluster_gas)
                a_bad = pop.badall(atom_all, H_fix)
                a_same_str = pop.is_same_structure(allatom_copy, atom_all)
                #a_same_str = 0
                if a_bad == 0 and a_same_str == 0:
                    x += 0
                    allatom = atom_all
                    cluster_gas = cluster_gas_new
                    with open('step-all', 'a') as fp:
                        fp.write(
                            '------%s-move-cluster-%s is finishing-------\n' % (times, i))
                    print('------%s-move-cluster-%s is finishing-------' %
                          (times, i))
                    break
                else:
                    if i < 1999:
                        allatom = allatom_copy
                        cluster_gas = cluster_gas_copy
                    else:
                        allatom = allatom_copy
                        cluster_gas = cluster_gas_copy
                        print('move-cluster-bad')
                        x += 1
                        break


# ###------------------get gas-----------------
        else:
            number_gas = len(gas_atom_ini)//2
            number2 = uniform(0, 1)
            #number2 = 1
            print(number_gas)
            number_cluster = len(cluster_atom)
            probality = 0.5
            if number_gas < 1:
                number2 = 0
            # elif number_gas >= number_cluster:
            #    number2 = 1

        # ###------------------get gas add-----------------
            if number2 < probality:
                for i in range(2000):
                    t = 1000 + times * 5
                    allatom_copy = copy.deepcopy(allatom)
                    cluster_gas_copy = copy.deepcopy(cluster_gas)
                    atom_all, cluster_gas_new = pop.add_gas(
                        support_atom, cluster_atom, t, allatom, cluster_gas)
                    a_bad = pop.badall(atom_all, H_fix)
                    if a_bad == 0:
                        x += 0
                        allatom = atom_all
                        cluster_gas = cluster_gas_new
                        with open('step-all', 'a') as fp:
                            fp.write(
                                '------%s-add-gas-%s is finishing-------\n' % (times, i))
                        print('------%s-add-gas-%s is finishing-------' %
                              (times, i))
                        break
                    else:
                        if i < 1999:
                            allatom = allatom_copy
                            cluster_gas = cluster_gas_copy
                        else:
                            allatom = allatom_copy
                            cluster_gas = cluster_gas_copy
                            print('add-gas-bad')
                            x += 1
                            break

        # ###------------------get gas del-----------------
            else:
                for i in range(2000):
                    allatom_copy = copy.deepcopy(allatom)
                    cluster_gas_copy = copy.deepcopy(cluster_gas)
                    atom_all, cluster_gas_post = pop.del_gas(
                        gas_atom_ini, allatom, cluster_gas)
                    a_bad = pop.badall(atom_all, H_fix)
                    if a_bad == 0:
                        x += 0
                        allatom = atom_all
                        cluster_gas = cluster_gas_post
                        with open('step-all', 'a') as fp:
                            fp.write(
                                '------%s-del gas-%s is finishing-------\n' % (times, i))
                        print('------%s-del gas-%s is finishing-------' %
                              (times, i))
                        break
                    else:
                        if i < 1999:
                            allatom = allatom_copy
                            cluster_gas = cluster_gas_copy
                        else:
                            allatom = allatom_copy
                            cluster_gas = cluster_gas_copy
                            print('del-gas-bad')
                            x += 1
                            break



        a = pop.badall(allatom, H_fix)

# #----------------------------------eann method----------------
        if x == 0 and a == 0:
            for i in range(2):
                support_atom_post, Height_post = pop.support_atom(
                    allatom, H_fix)
                cluster_atom_post = pop.cluster_atom(allatom, Height_post)
                gas_add = pop.gas_atom(
                    allatom, support_atom_post, cluster_atom_post)
                support_atom, Height = pop.support_atom(allatom, H_fix)
                cluster_gas, gas_no_choose = pop.cluster_adsorp(allatom, H_fix)
                if len(gas_no_choose) > 0:
                    allatom_copy = copy.deepcopy(allatom)
                    for pot in gas_no_choose:
                        for pot_gas in pot:
                            del allatom_copy[pot_gas]
                    allatom = allatom_copy
                atom_new = copy.deepcopy(allatom)
                atom_new = dict(
                    sorted(atom_new.items(), key=lambda d: d[0], reverse=False))
                for pot_fix in atom_new.keys():
                    if atom_new[pot_fix][2] > H_fix:
                        for i in range(3):
                            atom_new[pot_fix].append('T')
                    else:
                        for i in range(3):
                            atom_new[pot_fix].append('F')

                write_POSCAR(atom_new, 'POSCAR', 'POSCAR-1', H_fix)

                os.system('mv POSCAR-1 POSCAR')
                cell1 = ase.io.vasp.read_vasp("POSCAR")
                atoms = Atoms(cell1)
                traj = Trajectory('allposcar.traj', 'a')
                traj.write(atoms)

                a_bad_str = 0
                #cell1 = ase.io.vasp.read_vasp("POSCAR")
                #atoms = Atoms(cell1)
            # ----------------------eann --------------------------------
                atoms.calc = calc
                dyn = LBFGS(atoms, trajectory='atom2.traj')
                dyn.run(fmax=0.1, steps=100)
                traj = Trajectory('atom2.traj')
                os.system('mv atom2.traj opt/atom2-%s.traj' % times)
                atoms = traj[-1]
                ase.io.write('best-POSCAR', atoms, format='vasp', vasp5='True')
                if a_bad_str == 0:
                    #x_energy_bad = badstr_energy(Energy_all)
                    x_energy_bad = 0
                    allatom_gas = POSCAR_To_atom1(path, 'best-POSCAR')
                    support_atom_post, Height_post = pop.support_atom(
                        allatom_gas, H_fix)
                    cluster_atom_post = pop.cluster_atom(
                        allatom_gas, Height_post)
                    gas_post = pop.gas_atom(
                        allatom_gas, support_atom_post, cluster_atom_post)
                    ration_final = pop.cluster_ration(
                        cluster_atom_post, Height_post)
                    a_bad = pop.badall(allatom_gas, H_fix)
                    number_support = len(support_atom_post)
                    #number_CO2, allatom_no_CO2 = pop.CO2_bave(allatom_gas)
                    atoms.calc = calc
                    energy = atoms.get_potential_energy()
                    force_bad = atoms.get_forces()
                    ration_final = 0
                    number_CO2 = 0
                #    if number_CO2 == 0 :
                    #Energy = pop.energy_correction(num_atom_c, num_atom_O, 0, len(cluster_atom_post), energy)
                    #traj = Trajectory('allstr.traj', 'a')
                    # traj.write(atoms)
                    if ration_final == 0 and number_CO2 == 0 and a_bad == 0 and x_energy_bad == 0 and np.max(force_bad) <= 0.2:
                        #    Energy = pop.energy_correction(num_atom_c, num_atom_O, 0, len(cluster_atom_post), energy)
                        #    traj = Trajectory('allstr.traj', 'a')
                        #    traj.write(atoms)
                        Energy = pop.energy_correction(
                            num_atom_c, num_atom_O, 0, len(cluster_atom_post), energy)
                        traj = Trajectory('allstr.traj', 'a')
                        traj.write(atoms)
                        a = pop.energy_estimate(Energy, ini_energy)
                        if a == 0 or abs(Energy-best_energy) <= 0.5:
                            os.system('cp best-POSCAR bestPOSCAR')
                            os.system('mv best-POSCAR POSCAR')
                            traj1 = Trajectory('meta.traj', 'a')
                            traj1.write(atoms)
                            ini_energy = Energy
                            stop_time = cycle_time
                            print('add metastr')
                            if a == 0:
                                best_energy = Energy
                            break
                        else:
                            os.system('cp POSCAR-old POSCAR')
                            #os.system('rm best-POSCAR')
                            break
                    elif number_CO2 != 0:
                        if i <= 2:
                            # allatom = allatom_no_CO2
                            continue
                        else:
                            os.system('cp POSCAR-old POSCAR')
                            break
                    else:
                        os.system('cp POSCAR-old POSCAR')
                        break
                else:
                    os.system('cp POSCAR-old POSCAR')
                    break
        else:
            os.system('cp POSCAR-old POSCAR')
        if abs(stop_time - cycle_time) >= stop_criterion:
            print('finish')
            break
