import h5py
import numpy as np
import os
import glob
import pyvista as pv
import pandas as pd


def dump_h5_files(DFN, porosity):
    """Write permeability, porosity, and material ID values to HDF5 files for PFLOTRAN.

    Parameters
    ----------
    DFN : object
        DFN class instance; hydraulic properties must be set before calling.
    porosity : list
        Two-element list [matrix_porosity, fracture_porosity].

    Returns
    -------
    None
    """
    print('*' * 80)
    print("--> Dumping h5 file")

    fp_perm = 'permeability.h5'
    print(f'\n--> Opening HDF5 File {fp_perm}')
    with h5py.File(fp_perm, mode='w') as h5file:
        iarray = np.arange(1, DFN.num_nodes + 1)
        h5file.create_dataset('Cell Ids', data=iarray)
        print('--> Creating permeability array')
        for i in range(DFN.num_nodes):
            DFN.perm_cell[i] = DFN.perm[DFN.material_ids[i] - 1]
        h5file.create_dataset('Permeability', data=DFN.perm_cell)

    fp_porosity = 'porosity.h5'
    print(f'\n--> Opening HDF5 File {fp_porosity}')
    with h5py.File(fp_porosity, mode='w') as h5file:
        iarray = np.arange(1, DFN.num_nodes + 1)
        h5file.create_dataset('Cell Ids', data=iarray)
        porosity_cell = np.zeros_like(DFN.perm_cell)
        for i in range(DFN.num_nodes):
            porosity_cell[i] = porosity[DFN.material_ids[i] - 1]
        h5file.create_dataset('Porosity', data=porosity_cell)

    fp_matid = 'materials.h5'
    print(f'\n--> Opening HDF5 File {fp_matid}')
    with h5py.File(fp_matid, mode='w') as h5file:
        h5file.create_dataset('Materials/Cell Ids', data=iarray)
        h5file.create_dataset('Materials/Material Ids', data=DFN.material_ids)

    print("--> Done writing h5 files")
    print('*' * 80)
    print()


def get_uge_information(uge_filename='full_mesh.uge'):
    """Parse node coordinates and element volumes from a PFLOTRAN UGE file."""
    with open(uge_filename, 'r') as fuge:
        header = fuge.readline()
        num_cells = int(header.split()[-1])
        x   = np.zeros(num_cells)
        y   = np.zeros(num_cells)
        z   = np.zeros(num_cells)
        vol = np.zeros(num_cells)
        for i in range(num_cells):
            line = fuge.readline().split()
            x[i]   = float(line[1])
            y[i]   = float(line[2])
            z[i]   = float(line[3])
            vol[i] = float(line[4])
    return x, y, z, vol


def convert_uge_to_graph(uge_filename='full_mesh.uge'):
    """Build a NetworkX graph from UGE connectivity."""
    import networkx as nx
    G = nx.Graph()
    with open(uge_filename, 'r') as fuge:
        header = fuge.readline()
        num_cells = int(header.split()[-1])
        for _ in range(num_cells):
            line = fuge.readline().split()
            G.add_node(int(line[0]),
                       x=float(line[1]), y=float(line[2]),
                       z=float(line[3]), vol=float(line[4]))
        header = fuge.readline()
        num_conn = int(header.split()[-1])
        for _ in range(num_conn):
            line = fuge.readline().split()
            G.add_edge(int(line[0]), int(line[1]),
                       x=float(line[2]), y=float(line[3]),
                       z=float(line[4]), area=float(line[5]))
    return G


def save_results(index):
    """Save parsed VTK point data to CSV files (legacy path)."""
    output_dir = f"run_data_x{index:02d}"
    os.makedirs(output_dir, exist_ok=True)
    for vtk_file in glob.glob("parsed_vtk/*vtk"):
        print(f"Processing {vtk_file}")
        mesh = pv.read(vtk_file)
        if mesh.point_data.keys():
            df = pd.concat(
                [pd.DataFrame(mesh.points, columns=["X", "Y", "Z"]),
                 pd.DataFrame({n: mesh.point_data[n] for n in mesh.point_data.keys()})],
                axis=1)
            base_name = os.path.splitext(os.path.basename(vtk_file))[0]
            df.to_csv(os.path.join(output_dir, f"{base_name}_point_data.csv"), index=False)
        else:
            print("Warning: no point data arrays found.")



def save_results_h5(index, x, y, path, params=None):
    """Write pressure and permeability fields to an HDF5 file.
 
    The HDF5 file contains:
      - Dataset  "grid_points"   : (N,2) float64 array of (x, y) node coords
      - Dataset  "pressure"      : (N,)  float64 liquid pressure [Pa]
      - Dataset  "permeability"  : (N,)  float64 cell permeability [m^2]
      - Attributes on root group : contents of `params` dict (scalar floats)
 
    Parameters
    ----------
    index : int
        Simulation index used to construct file names.
    x, y  : array-like
        Node x- and y-coordinates from get_uge_information().
    path  : str
        Directory where the output HDF5 file is written.
    params : dict, optional
        Physical / sweep parameters to store as HDF5 root attributes, e.g.::
 
            {
                'k_f':          1e-8,   # fracture permeability  [m^2]
                'k_m':          1e-16,  # matrix permeability    [m^2]
                'k_ratio':      1e8,    # k_f / k_m              [-]
                'p32_target':   1.0,    # target fracture intensity [m^-1]
                'p21_realized': 0.87,   # realized 2-D intensity  [m^-1]
            }
 
        Any additional key/value pairs are written verbatim.
    """
    os.makedirs(path, exist_ok=True)
    h5_path      = os.path.join(path, f"pressure_x{index:02d}.h5")
    pflotran_h5  = f"pflotran_lhs_{index:02d}_pressure.h5"
 
    # PFLOTRAN HDF5 groups are named "0 Time 0.00000E+00 y" (IC) and
    # "1 Time 1.00000E+03 y" (steady state).  Sort keys and take the last one.
    with h5py.File(pflotran_h5, "r") as src:
        final_group = sorted(src.keys())[-1]
        pressure    = src[final_group]["Liquid Pressure [Pa]"][:]
        permeability = src[final_group]["Permeability [m^2]"][:]
 
    with h5py.File(h5_path, "w") as h5f:
        # ── Physical-parameter header ────────────────────────────────────────
        if params:
            for key, val in params.items():
                h5f.attrs[key] = float(val)
            print(f"--> H5 header attributes written: {list(params.keys())}")
 
        # ── Grid and field data ──────────────────────────────────────────────
        h5f.create_dataset("grid_points",  data=np.column_stack([x, y]),
                           dtype=np.float64)
        h5f.create_dataset("pressure",     data=pressure,    dtype=np.float64)
        h5f.create_dataset("permeability", data=permeability, dtype=np.float64)
 
    print(f"--> All results written to {h5_path}")
