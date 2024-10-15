# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 15:52:34 2021

@author: DELL
"""

import ase.io.vasp
from ase import Atoms
from ase.calculators.eann import EANN
from ase.optimize.minimahopping import MinimaHopping
from ase.optimize.minimahopping import MHPlot
from ase.optimize import BFGS, FIRE, LBFGS
import os
#from ase.neb import NEB
from ase.autoneb import AutoNEB 
from ase.optimize import MDMin, BFGS
from ase.io import read
from ase.dyneb import DyNEB
from ase.io.trajectory import Trajectory
from ase.io import read, write
import sys

atomtype = ['Cu','O','C']
nn = 'EANN_PES_DOUBLE.pt'
pos1 = sys.argv[1]
pos2 = sys.argv[2]
os.system('rm prefix*.traj;rm -r AutoNEB_iter')
#geom_ini = ase.io.vasp.read_vasp("POSCAR1")
#geom_fin = ase.io.vasp.read_vasp("POSCAR2-2")
geom_ini = ase.io.vasp.read_vasp("%s" %pos1)
geom_fin = ase.io.vasp.read_vasp("%s" %pos2)
geom_ini.calc = EANN(atomtype=atomtype,period=[1,1,1],nn=nn)
geom_ini.get_potential_energy()
geom_fin.calc = EANN(atomtype=atomtype,period=[1,1,1],nn=nn)
geom_fin.get_potential_energy()
write('prefix000.traj', geom_ini)
write('prefix001.traj', geom_fin)
#traj = Trajectory('prefix000.traj', 'w', geom_ini)
#traj = Trajectory('prefix001.traj', 'w', geom_fin)
#geom_ini = read('atom1.traj')
#geom_fin = read('atom2.traj')
#nimg = 10
calc = EANN(atomtype=atomtype,period=[1,1,1],nn=nn)
#images = [geom_ini.copy() for i in range(nimg + 1)] + [geom_fin]
#images = ase.io.read('A2B.traj@-12:')
#images = [initial]
#for i in range(nimg):
#    image = initial.copy()
#    image.calc = EANN(atomtype=atomtype,period=[1,1,0],nn=nn)
#    images.append(image)
#images.append(final)
#for img in images:
#    img.calc = EANN(atomtype=atomtype,period=[1,1,0],nn=nn)
#    print(img)
#geom_ini.calc = calc
#geom_fin.calc = calc
#-----------------------auto neb--------------------------------
band = AutoNEB(attach_calculators=calc, # k=0.05,
               prefix='prefix',n_simul=1,n_max=6,fmax=0.1, climb=True, smooth_curve=True)
               #remove_rotation_and_translation=True, smooth_curve=True,)
band.run()
#band.optimizer.run(fmax=0.1)
#band = DyNEB(images, fmax=0.1, dynamic_relaxation=True)
#band = NEB(images,climb=True)
#band.interpolate(method='idpp')
# Interpolate linearly the potisions of the three middle images:
#optimizer = MDMin(band, trajectory='neb.traj')
#optimizer = LBFGS(band, trajectory='A2B2.traj')
#optimizer.run(fmax=0.2)
