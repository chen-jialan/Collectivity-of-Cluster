from ase.ga.data import DataConnection
import ase.io.vasp
from ase.io import write
from ase.db import connect
import copy
import os
import numpy as np

if os.path.exists('./bulk.db'):
    os.remove('./bulk.db')
if os.path.exists('./pos'):
    os.system('rm -r ./pos')
os.mkdir('./pos')
da = connect('gadb.db')
energy = [da.get(id=i).energy for i in range(1,len(da)+1) if 'energy' in da.get(id=i)]
index_atom = [i for i in range(1,len(da)+1) if 'energy' in da.get(id=i)]
               #da.get(id=i).energy >= -1320.653819639333]
energy_copy = copy.deepcopy(energy)
energy_copy.sort()
index_copy = [index_atom[energy.index(e)] for e in energy_copy]
#print(energy_copy,index_copy)
with open ('./energy_ga_all',"w") as f1:
    for i  in range(len(energy_copy)):
        f1.write('%s %s\n' %(energy_copy[i],index_copy[i]))
f1.close()

#for i in range(100)
atom_all = []
a = 0
f = open('./energy',"w")
energy_write = min(energy_copy)
for i in range(len(energy_copy)):
    if abs(energy_copy[i]-min(energy_copy)) <= 1:
        atoms = da.get(id=index_copy[i]).toatoms()
        if len(atom_all) == 0:
            atom_all.append(atoms)
            ase.io.write('POSCAR-final', atoms, format='vasp', vasp5='True')
            os.system('mv POSCAR-final pos/POSCAR-%s' %a)
            f.write('%s %s\n' %(a,energy_copy[i]))
            a += 1
        else:
            j = 1
            a2_add_position = atoms.get_positions()
            a2_std = np.std(a2_add_position)
            for atom_e in atom_all:
                a1_position = atom_e.get_positions()
                a1_std = np.std(a1_position)
                #print(np.std(a1_position-a2_add_position))
                print(a1_std, a2_std)
                if abs(a1_std-a2_std) <=2e-4:
                    j =0
                    break
            if j == 1 and abs(energy_copy[i]-energy_write) >= 0.01:
                print('add structure')
                ase.io.write('POSCAR-final', atoms, format='vasp', vasp5='True')
                os.system('mv POSCAR-final pos/POSCAR-%s' %a) 
                atom_all.append(atoms)
                f.write('%s %s\n' %(a,energy_copy[i]))
                energy_write = energy_copy[i]
                a += 1
        # print(atoms.get_potential_energy())
os.system('mv energy pos/.')
f.close()
db = connect('bulk.db')
for atoms in atom_all:
    db.write(atoms, relaxed=True)
        #db = connect('bulk.db')
        #db.write(atoms, relaxed=True)
