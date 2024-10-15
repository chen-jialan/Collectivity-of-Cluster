# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 15:52:34 2021

@author: DELL
"""

import ase.io.vasp
from ase import Atoms
from ase.io.trajectory import Trajectory
import numpy as np
import os

#cell1 = ase.io.vasp.read_vasp("POSCAR_3")
#atoms = Atoms(cell1)
traj = Trajectory('meta.traj')
#traj.write(atoms)
for i in range(len(traj)):
    ase.io.write('POSCAR-meta', traj[i], format='vasp', vasp5='True') 
    os.system('cp POSCAR-meta pos/POSCAR-%s' %i)
