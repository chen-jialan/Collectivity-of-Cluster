from ase.io.trajectory import Trajectory
import os
import ase.io.vasp

path = os.getcwd()
file_name = 'nebtraj'
file_traj = 'A2B.traj'
file_all = [file_e for file_e in os.listdir(path) if 'pre' in file_e]
number = [str(i).rjust(3,'0') for i in range(len(file_all))]


if os.path.exists(file_name):
    os.system('rm -r %s' %file_name)
if os.path.exists(file_name):
    os.remove(file_traj)
os.mkdir('./%s' %file_name)


for i in range(len(number)):
    file_e  = 'prefix%s.traj'  %number[i]
    print(file_e) 
    traj = Trajectory('%s' %file_e)[-1]
    ase.io.write('./%s/POSCAR-%s' %(file_name,i), traj, format='vasp', vasp5='True') 
    traj_write = Trajectory(file_traj, 'a') 
    traj_write.write(traj) 
