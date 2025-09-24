""" DFM Driver """

from pydfnworks import *
import os, subprocess, sys
import numpy as np
from helper_dump_h5 import dump_h5_files, convert_uge_to_graph, save_results, save_results_h5, get_uge_information
from write_pflotran_card_2d import write_pflotran_card_pressure

import pandas as pd
import shutil 

p32 = 1

domain_index = int(sys.argv[1])
lhs_index = domain_index - 1

src_path = os.getcwd()
jobname = f"{src_path}/pressure_x{domain_index:02d}"
dfnFlow_file = f"{src_path}/pflotran_files/pflotran_lhs_{domain_index:02d}_pressure.in"
write_pflotran_card_pressure(dfnFlow_file,2e6)

DFN = DFNWORKS(jobname, dfnFlow_file=dfnFlow_file, ncpu=1)

DFN.params['domainSize']['value'] = [10, 10, 10]
DFN.params['domainSizeIncrease']['value'] = [2, 2, 0]
DFN.params['h']['value'] = 0.1
DFN.params['tripleIntersections']['value'] = True
DFN.params['stopCondition']['value'] = 1
DFN.params['ignoreBoundaryFaces']['value'] = True
DFN.params['keepOnlyLargestCluster']['value'] = False
DFN.params['seed']['value'] = domain_index*10
DFN.params['rFram']['value'] = True

DFN.add_user_fract_from_file(
    shape="poly",
    filename=f'{src_path}/domain.dat',
    permeability=1 * [1e-12],  #list or array of nPolygons perms
    nPolygons=1)

DFN.add_fracture_family(
    shape="rect",
    distribution="tpl",
    alpha=2.4,
    min_radius=3.0,
    max_radius=15.0,
    kappa=20.0,
    theta=90.0,
    phi=0.0,
    p32=p32,
    hy_variable='aperture',
    hy_function='correlated',
    hy_params={
        "alpha": 10**-5,
        "beta": 0.5
    })

DFN.add_fracture_family(
    shape="rect",
    distribution="tpl",
    alpha=2.4,
    min_radius=3.0,
    max_radius=15.0,
    kappa=20.0,
    theta=90.0,
    phi=90.0,
    #aspect=2,
    p32=p32,
    hy_variable='aperture',
    hy_function='correlated',
    hy_params={
        "alpha": 10**-5,
        "beta": 0.5
    })

DFN.make_working_directory(delete=True)

DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()
DFN.num_frac = 1
DFN.mesh_network(uniform_mesh=True, strict=False)

cmd = f'lagrit < {src_path}/process_mesh.lgi'
DFN.call_executable(cmd)
DFN.aperture = 10 * np.ones(DFN.num_frac)

DFN.lagrit2pflotran()
DFN.zone2ex(zone_file='boundary_left.zone', face='west')
DFN.zone2ex(zone_file='boundary_right.zone', face='east')
DFN.uge_file = 'full_mesh.uge'

DFN.material_ids = np.genfromtxt('materialid.dat', skip_header=3).astype(int)

## Settting perm / porosity parameters
DFN.perm[0] = 1e-16
DFN.perm[1] = 1e-8
porosity = [0.01, 0.5]

dump_h5_files(DFN, porosity)

matrix_cells = np.where(DFN.material_ids == 1)[0]
fracture_cells = [np.where(DFN.material_ids == 2)[0]]

with open('matrix.txt', 'w') as fp:
    for i in matrix_cells:
        fp.write(f"{i+1}\n")

with open('fracture.txt', 'w') as fp:
    for i in fracture_cells[0]:
        fp.write(f"{i+1}\n")


output_dir = f'run_data_x{domain_index:02d}' 
os.makedirs(output_dir)

DFN.ncpu = 4
DFN.pflotran()
DFN.parse_pflotran_vtk_python()
x,y,z,volume = get_uge_information()
save_results_h5(domain_index, x,y, src_path + os.sep + 'h5_files')  