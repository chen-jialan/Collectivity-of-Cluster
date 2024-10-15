from ase.ga.data import DataConnection
from ase.io import write
from ase.io.trajectory import Trajectory

da = DataConnection('gadb.db')
write('all_candidates.traj', da.get_all_relaxed_candidates())
traj = Trajectory('all_candidates.traj')
for i in range(len(traj)):
    print(traj[i].get_potential_energy()) 
