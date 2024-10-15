import matplotlib.pyplot as plt
from ase.neb import NEBTools
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.calculators.eann import EANN
import ase.io.vasp
import os

if os.path.exists('neb_traj'):
    os.system('rm -r neb_traj')
os.mkdir('./neb_traj')

images = read('A2B2.traj@-6:')
os.system('cp A2B2.traj ./neb_traj/.')
#atoms = Trajectory('A2B.traj')[-1]
nebtools = NEBTools(images)
atomtype = ['Cu','O','C']
nn = 'EANN_PES_DOUBLE.pt'
#image.calc = EANN(atomtype=atomtype,period=[1,1,0],nn=nn)
i = 0
for imag in images:
    imag.calc = EANN(atomtype=atomtype,period=[1,1,1],nn=nn)
    ase.io.write('./neb_traj/POSCAR-%s' %i, imag, format='vasp', vasp5='True') 
    i += 1

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()
# Get the barrier without any interpolation between highest images.
Ef, dE = nebtools.get_barrier(fit=False)

# Get the actual maximum force at this point in the simulation.
max_force = nebtools.get_fmax()

# Create a figure like that coming from ASE-GUI.
fig = nebtools.plot_band()
fig.savefig('diffusion-barrier.png')

# Create a figure with custom parameters.
fig = plt.figure(figsize=(5.5, 4.0))
ax = fig.add_axes((0.15, 0.15, 0.8, 0.75))
nebtools.plot_band(ax)
fig.savefig('diffusion-barrier.png')
