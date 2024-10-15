# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 15:52:34 2021

@author: DELL
"""

import ase.io.vasp
from ase import Atoms
from ase.calculators.eann_version1 import EANN
from ase.optimize.minimahopping import MinimaHopping
from ase.optimize.minimahopping import MHPlot
from ase.optimize import LBFGS
from ase.neb import NEB
#from ase.optimize import MDMin, BFGS
from ase.io import read
from ase.dyneb import DyNEB
import torch


geom_ini = ase.io.vasp.read_vasp("POSCAR1")
geom_fin = ase.io.vasp.read_vasp("POSCAR2")
#geom_ini = read('atom1.traj')
#geom_fin = read('atom2.traj')
nimg = 4
atomtype = ['Cu', 'O', 'C']
nn = 'EANN_PES_DOUBLE.pt'
pes = torch.jit.load(nn)
pes.to('cpu').to(torch.double)
pes.eval()
images = [geom_ini.copy() for i in range(nimg + 1)] + [geom_fin]
#images = ase.io.read('A2B.traj@-12:')
#images = [initial]
for img in images:
    img.calc = EANN(atomtype=atomtype, period=[1, 1, 0], pes=pes)
#    print(img)
band = DyNEB(images, fmax=0.1, dynamic_relaxation=True, scale_fmax=2.)
#band = NEB(images,climb=True)
band.interpolate(method='idpp')
# Interpolate linearly the potisions of the three middle images:
#optimizer = MDMin(band, trajectory='neb.traj')
optimizer = LBFGS(band, trajectory='A2B2.traj')
optimizer.run(fmax=0.1)
