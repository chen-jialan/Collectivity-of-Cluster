# -*- coding: utf-8 -*-
"""
NEB or DyNEB Runner Script with Modular Design

Author: Jia-Lan Chen
Date: 2025-08-02
python dyneb.py --poscar1 POSCAR1 --poscar2 POSCAR2 \
  --pes EANN_PES_DOUBLE.pt --atomtype Cu O C \
  --n_images 5 --traj myneb.traj --dyneb
"""

import argparse
import ase.io
from ase.neb import NEB
from ase.dyneb import DyNEB
from ase.optimize import LBFGS
from ase.io import read
from ase.calculators.eann_version1 import EANN
import torch

def load_calculator(atomtype, pes_path, period):
    pes = torch.jit.load(pes_path)
    pes.to('cpu').to(torch.double)
    pes.eval()
    return lambda atoms: EANN(atomtype=atomtype, period=period, pes=pes)

def build_images(geom_ini, geom_fin, n_images, calculator_func):
    images = [geom_ini.copy() for _ in range(n_images + 1)] + [geom_fin]
    for img in images:
        img.calc = calculator_func(img)
    return images

def run_dyneb(
    poscar1_path='POSCAR1',
    poscar2_path='POSCAR2',
    pes_path='EANN_PES_DOUBLE.pt',
    atomtype=None,
    period=None,
    n_images=4,
    fmax=0.1,
    traj_file='neb_output.traj',
    use_dyneb=True
):
    if atomtype is None:
        atomtype = ['Cu', 'O', 'C']
    if period is None:
        period = [1, 1, 0]

    geom_ini = ase.io.read(poscar1_path)
    geom_fin = ase.io.read(poscar2_path)

    calc_func = load_calculator(atomtype, pes_path, period)
    images = build_images(geom_ini, geom_fin, n_images, calc_func)

    neb_class = DyNEB if use_dyneb else NEB
    band = neb_class(images, fmax=fmax, dynamic_relaxation=True, scale_fmax=2.) if use_dyneb else neb_class(images, climb=True)
    band.interpolate(method='idpp')

    optimizer = LBFGS(band, trajectory=traj_file)
    optimizer.run(fmax=fmax)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Modular NEB or DyNEB Runner")

    parser.add_argument('--poscar1', type=str, default='POSCAR1', help='Initial structure file (VASP format)')
    parser.add_argument('--poscar2', type=str, default='POSCAR2', help='Final structure file (VASP format)')
    parser.add_argument('--pes', type=str, default='EANN_PES_DOUBLE.pt', help='Path to the PES model')
    parser.add_argument('--atomtype', nargs='+', default=['Cu', 'O', 'C'], help='List of atom types used in the system')
    parser.add_argument('--period', nargs=3, type=int, default=[1, 1, 0], help='Periodicity in x, y, z directions')
    parser.add_argument('--n_images', type=int, default=4, help='Number of intermediate images')
    parser.add_argument('--fmax', type=float, default=0.1, help='Maximum force criteria for convergence')
    parser.add_argument('--traj', type=str, default='neb_output.traj', help='Trajectory file to save')
    parser.add_argument('--dyneb', action='store_true', help='Use DyNEB instead of regular NEB')

    args = parser.parse_args()

    run_dyneb(
        poscar1_path=args.poscar1,
        poscar2_path=args.poscar2,
        pes_path=args.pes,
        atomtype=args.atomtype,
        period=args.period,
        n_images=args.n_images,
        fmax=args.fmax,
        traj_file=args.traj,
        use_dyneb=args.dyneb
    )
