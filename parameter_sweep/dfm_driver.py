"""DFM Driver — reads per-sample parameters from lhs_params.csv."""

from pydfnworks import *
import os, sys
import numpy as np
import pandas as pd
from helper_dump_h5 import dump_h5_files, save_results_h5, get_uge_information
from write_pflotran_card_2d import write_pflotran_card_pressure

# ── Sample index from command line ────────────────────────────────────────────
domain_index = int(sys.argv[1])

# ── Load sweep parameters ─────────────────────────────────────────────────────
lhs_df      = pd.read_csv('lhs_params.csv', index_col='domain_index')
k_f         = float(lhs_df.loc[domain_index, 'k_f'])
k_m         = float(lhs_df.loc[domain_index, 'k_m'])
k_ratio     = float(lhs_df.loc[domain_index, 'k_ratio'])
p32_val     = float(lhs_df.loc[domain_index, 'p32'])
sample_type = str(lhs_df.loc[domain_index, 'sample_type'])

print(f"--> Sample {domain_index:04d} [{sample_type}]: "
      f"k_m={k_m:.2e}  k_ratio={k_ratio:.1e}  k_f={k_f:.2e}  p32={p32_val:.3f}")

# ── Paths ─────────────────────────────────────────────────────────────────────
src_path     = os.getcwd()
jobname      = f"{src_path}/pressure_x{domain_index:02d}"
dfnFlow_file = f"{src_path}/pflotran_files/pflotran_lhs_{domain_index:02d}_pressure.in"
write_pflotran_card_pressure(dfnFlow_file, 2e6)

# ── DFN setup ─────────────────────────────────────────────────────────────────
DFN = DFNWORKS(jobname, dfnFlow_file=dfnFlow_file, ncpu=1)

DFN.params['domainSize']['value']             = [10, 10, 10]
DFN.params['domainSizeIncrease']['value']     = [2, 2, 0]
DFN.params['h']['value']                      = 0.1
DFN.params['tripleIntersections']['value']    = True
DFN.params['stopCondition']['value']          = 1
DFN.params['ignoreBoundaryFaces']['value']    = True
DFN.params['keepOnlyLargestCluster']['value'] = False
DFN.params['seed']['value']                   = domain_index * 10
DFN.params['rFram']['value']                  = True

DFN.add_user_fract_from_file(
    shape="poly",
    filename=f'{src_path}/domain.dat',
    permeability=[1e-12],
    nPolygons=1)

for phi in [0.0, 90.0]:
    DFN.add_fracture_family(
        shape="rect",
        distribution="tpl",
        alpha=2.4,
        min_radius=3.0,
        max_radius=15.0,
        kappa=20.0,
        theta=90.0,
        phi=phi,
        p32=p32_val,
        hy_variable='aperture',
        hy_function='correlated',
        hy_params={"alpha": 10**-5, "beta": 0.5})

DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()

# ── Realized p21 from DFN object ──────────────────────────────────────────────
p21_realized = float(DFN.p21[0])
print(f"--> Realized p21 = {p21_realized:.4f} m^-1  (target p32 = {p32_val:.4f} m^-1)")

DFN.num_frac = 1
DFN.mesh_network(uniform_mesh=True, strict=False)

cmd = f'lagrit < {src_path}/process_mesh.lgi'
DFN.call_executable(cmd)
DFN.aperture = 10 * np.ones(DFN.num_frac)

DFN.lagrit2pflotran(dim  = 2 )
DFN.zone2ex(zone_file='boundary_left.zone',  face='west')
DFN.zone2ex(zone_file='boundary_right.zone', face='east')
DFN.uge_file = 'full_mesh.uge'

DFN.material_ids = np.genfromtxt('materialid.dat', skip_header=3).astype(int)

# ── Hydraulic properties ──────────────────────────────────────────────────────
DFN.perm[0] = k_m
DFN.perm[1] = k_f
porosity = [0.01, 0.5]
dump_h5_files(DFN, porosity)

matrix_cells   = np.where(DFN.material_ids == 1)[0]
fracture_cells = np.where(DFN.material_ids == 2)[0]

with open('matrix.txt', 'w') as fp:
    for i in matrix_cells:
        fp.write(f"{i+1}\n")
with open('fracture.txt', 'w') as fp:
    for i in fracture_cells:
        fp.write(f"{i+1}\n")

output_dir = f'run_data_x{domain_index:02d}'
os.makedirs(output_dir, exist_ok=True)

# ── PFLOTRAN solve ────────────────────────────────────────────────────────────
DFN.ncpu = 4
DFN.pflotran()
# DFN.parse_pflotran_h5()

x, y, z, volume = get_uge_information()

# ── Save H5 with full parameter header ───────────────────────────────────────
params = {
    'k_f':          k_f,
    'k_m':          k_m,
    'k_ratio':      k_ratio,
    'p32_target':   p32_val,
    'p21_realized': p21_realized,
    'sample_type':  0 if sample_type == 'lhs' else 1,   # 0=lhs, 1=diagnostic
}
save_results_h5(domain_index, x, y, src_path + os.sep + 'h5_files', params=params)
