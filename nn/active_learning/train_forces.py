# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 00:49:34 2022

@author: jlchen
"""

import os
from NN_train import NSDS_algorithm, NN_predicted
import ase.io.vasp
from ase.io.trajectory import Trajectory
from ase.db import connect
import ase.io
import numpy as np

filename = 'allstr.traj'
atomtype1 = ['Cu', "C", "O"]
atomtype2 = atomtype1
nn1 = 'EANN_PES_DOUBLE.pt'
nn2 = 'EANN_PES_DOUBLE2.pt'
atom_nn_traj = NN_predicted(filename=filename, filetype='traj',
                            atomtype1=atomtype1, atomtype2=atomtype2, nn1=nn1, nn2=nn2)

atoms_all = atom_nn_traj.force_convergence(force_min=0.2, force_max=20)
ase.io.write('nn_train.traj', atoms_all)
