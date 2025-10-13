""" DFM Driver """

from pydfnworks import *
import os, subprocess, sys
import numpy as np
from helper_dump_h5 import dump_h5_files, convert_uge_to_graph, save_results, save_results_h5, get_uge_information
from write_pflotran_card_2d import write_pflotran_card_pressure

import pandas as pd
import shutil 
import random 

def write_domain(index):
    random.seed(index)
    y = random.uniform(-4, 4)
    filename = f'domains/domain_x{index}.dat'
    with open(filename, 'w') as fp:
        polygon_string = f"""nPolygons: 2
    4 {{-5,-5,0}} {{5,-5,0}} {{5, 5, 0}} {{-5, 5, 0}}
    4 {{-5,{y:.3f},-5}} {{-5,{y:.3f},5}} {{5,{y:.3f},5}} {{5,{y:.3f},-5}}
    """
        fp.write(polygon_string)

    return filename 
domain_index = int(sys.argv[1])
domain_filename = write_domain(domain_index)

src_path = os.getcwd()
jobname = f"{src_path}/single_fracture_pressure_x{domain_index}"
dfnFlow_file = f"{src_path}/pflotran_files/pflotran_lhs_{domain_index}_pressure.in"
write_pflotran_card_pressure(dfnFlow_file,2e6)

DFN = DFNWORKS(jobname, dfnFlow_file=dfnFlow_file, ncpu=1)

DFN.params['domainSize']['value'] = [10, 10, 10]
DFN.params['domainSizeIncrease']['value'] = [2, 2, 0]
DFN.params['h']['value'] = 0.5
DFN.params['tripleIntersections']['value'] = True
DFN.params['stopCondition']['value'] = 0
DFN.params['nPoly']['value'] = 2
DFN.params['ignoreBoundaryFaces']['value'] = True
DFN.params['keepOnlyLargestCluster']['value'] = False
DFN.params['seed']['value'] = domain_index*10
DFN.params['rFram']['value'] = True

DFN.add_user_fract_from_file(
    shape="poly",
    filename=f'{src_path}/{domain_filename}',
    permeability=2 * [1e-12],  #list or array of nPolygons perms
    nPolygons=2)

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