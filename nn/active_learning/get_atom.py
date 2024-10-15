# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 15:00:11 2022

@author: DELL
"""
import sys
from ase.db import connect
import ase.io.vasp
import os
import ase.io

if __name__ == '__main__':
   # cycle = int(sys.argv[1])
    if os.path.exists('./pos'):
        os.system('rm -r ./pos')
    os.mkdir('./pos')
    filename = 'nn_train.traj'
    #db = connect(filename)
    #print(len(db))
    atoms_all = ase.io.read(filename,index=":")
    for i in range(len(atoms_all)):
        atom = atoms_all[i]
        #atom = db.get(id=i+1).toatoms()
        ase.io.write('POSCAR', atom, format='vasp', vasp5='True' )
        os.system('mv POSCAR pos/POSCAR-%s' %(i))
        print(atom)
    
        
