# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 12:37:18 2025

@author: Administrator
"""

from gcmc import run
from ase.calculators.eann import EANN

if __name__ == "__main__":

    # --------------the type of atom--------------
    atomtype = ['Cu', 'Ce', 'O', 'C']
    # -----------------the device is cpu or gpu( cpu is default)---------------
    device = 'cpu'
    # -------------------------pbc([1,1,1] is default)---------------------
    period = [1, 1, 1]
    # ---------------nn file('EANN_PES_DOUBLE.pt' is default)----------------------------
    nn = 'EANN_PES_DOUBLE.pt'

    calc = EANN(device=device, nn=nn, atomtype=atomtype, period=period)
    # Reference molecular energies (in eV)
    potential_O2 = -1.071316E+01/2    # O2 energy (per atom)
    potential_CO = -15.78612          # CO energy
    potential_CO2 = -24.24603         # CO2 energy
    ratio_CO = 0.5                    # ratio pressure of CO
    # Simulation parameters
    temperature = 300      # Temperature (in K)
    pressure = 1e5         # Pressure (in Pa)
    cycle = 1000           # Total simulation steps
    # Early stopping criterion (stop if no change for 300 steps)
    stop_criterion = 300
    H_fix = 1.6           # Fixed height for support atoms (in Å)

    # Run the simulation
    run(calc=calc,                     # calculator
        potential_O2=potential_O2,     # O2 reference energy
        potential_CO=potential_CO,     # CO reference energy
        potential_CO2=potential_CO2,   # CO2 reference energy
        temperature=temperature,       # Simulation temperature
        pressure=pressure,            # Simulation pressure
        cycle=cycle,                  # Total simulation cycles
        stop_criterion=stop_criterion,  # Stopping condition
        H_fix=H_fix,                   # constraint
        ratio_CO=ratio_CO              # ratio pressure of CO
        )
