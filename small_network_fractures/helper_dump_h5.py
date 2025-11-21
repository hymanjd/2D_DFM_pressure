import h5py
import numpy as np 
import h5py

def dump_h5_files(DFN, porosity):
    """ Write permeability values to cell ids and permeability values to dfn_properties.h5 file for pflotran. 

    Parameters
    ----------
        self : object
            DFN Class

    Returns
    ---------
        None

    Notes
    ----------
        Hydraulic properties need to attached to the class prior to running this function. Use DFN.assign_hydraulic_properties() to do so. 
    """
    print('*' * 80)
    print("--> Dumping h5 file")
    fp_perm = 'permeability.h5'
    print(f'\n--> Opening HDF5 File {fp_perm}')
    with h5py.File(fp_perm, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        iarray = np.arange(1,DFN.num_nodes + 1)
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print('--> Creating permeability array')
        print('--> Note: This script assumes isotropic permeability')
        for i in range(DFN.num_nodes):
            DFN.perm_cell[i] = DFN.perm[DFN.material_ids[i] - 1]
        print('--> Writting Permeability')
        dataset_name = 'Permeability'
        h5dset = h5file.create_dataset(dataset_name, data=DFN.perm_cell)

    fp_porosity = 'porosity.h5'
    print(f'\n--> Opening HDF5 File {fp_porosity}')
    with h5py.File(fp_porosity, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        iarray = np.arange(1,DFN.num_nodes + 1)
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print('--> Creating porosity array')
        porosity_cell = np.zeros_like(DFN.perm_cell)
        for i in range(DFN.num_nodes):
            porosity_cell[i] = porosity[DFN.material_ids[i] - 1]
        print('--> Writting Porosity')
        dataset_name = 'Porosity'
        h5dset = h5file.create_dataset(dataset_name, data=porosity_cell)

    fp_matid = 'materials.h5'
    print(f'\n--> Opening HDF5 File {fp_matid}')
    with h5py.File(fp_matid, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        dataset_name = 'Materials/Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        dataset_name = 'Materials/Material Ids'
        h5dset = h5file.create_dataset(dataset_name, data=DFN.material_ids)


    print("--> Done writting h5 file")
    print('*' * 80)
    print()


def get_uge_information(uge_filename = 'full_mesh.uge'):
    with open(uge_filename, 'r') as fuge:
        header = fuge.readline()
        num_cells = int(header.split()[-1])
        x = np.zeros(num_cells)
        y = np.zeros(num_cells)
        z = np.zeros(num_cells)
        vol = np.zeros(num_cells)
        for i in range(num_cells):
            line = fuge.readline()
            line = line.split()
            # id = int(line[0])
            x[i] = float(line[1])
            y[i] = float(line[2])
            z[i] = float(line[3])
            vol[i] = float(line[4])
    return x, y, z, vol 

def convert_uge_to_graph(uge_filename = 'full_mesh.uge'):
    import networkx as nx
    G = nx.Graph()
    with open(uge_filename, 'r') as fuge:
        header = fuge.readline()
        num_cells = int(header.split()[-1])
        for _ in range(num_cells):
            line = fuge.readline()
            line = line.split()
            id = int(line[0])
            x = float(line[1])
            y = float(line[2])
            z = float(line[3])
            vol = float(line[4])
            G.add_node(id, x=x, y=y, z=z, vol=vol)
        header = fuge.readline()
        num_conn = int(header.split()[-1])
        for _ in range(num_conn):
            line = fuge.readline()
            line = line.split()
            id1 = int(line[0])
            id2 = int(line[1])
            x = float(line[2])
            y = float(line[3])
            z = float(line[4])  
            area = float(line[5])
            G.add_edge(id1, id2, x=x, y=y, z=z, area=area)

    return G 


import pyvista as pv
import pandas as pd
import os
import glob

def save_results(index):
    # Output directory for CSV files
    output_dir = f"run_data_x{index:02d}"
    os.makedirs(output_dir, exist_ok=True)
    vtk_files = glob.glob("parsed_vtk/*vtk")
    # Loop through each file and export cell data
    for vtk_file in vtk_files:
        print(f"Processing {vtk_file}")
        mesh = pv.read(vtk_file)

        # Print available point data arrays
        print("Available POINT data arrays:", mesh.point_data.keys())

        # Convert point data to DataFrame
        if mesh.point_data.keys():
            df = pd.DataFrame({name: mesh.point_data[name] for name in mesh.point_data.keys()})

            # Optional: include spatial coordinates
            coords = pd.DataFrame(mesh.points, columns=["X", "Y", "Z"])
            df = pd.concat([coords, df], axis=1)

            # Save CSV
            base_name = os.path.splitext(os.path.basename(vtk_file))[0]
            csv_path = os.path.join(output_dir, f"{base_name}_point_data.csv")
            df.to_csv(csv_path, index=False)
            print(f"Saved: {csv_path}")
        else:
            print("Warning. No point data arrays found.")


# def save_results_h5(index):
#     # Output directory for CSV files
#     output_dir = f"run_data_x{index:02d}"
#     os.makedirs(output_dir, exist_ok=True)
#     vtk_files = glob.glob("parsed_vtk/*vtk")
#     # Loop through each file and export cell data
#     for vtk_file in vtk_files:
#         print(f"Processing {vtk_file}")
#         mesh = pv.read(vtk_file)

#         # Get base filename
#         base_name = os.path.splitext(os.path.basename(vtk_file))[0]
#         h5_path = os.path.join(output_dir, f"{base_name}_point_data.h5")

#         # Open HDF5 file for writing
#         with h5py.File(h5_path, "w") as h5f:
#             # Save coordinates
#             h5f.create_dataset("points", data=mesh.points)

#             # Save point data arrays
#             for name in mesh.point_data.keys():
#                 h5f.create_dataset(name, data=mesh.point_data[name])

#         print(f"✔️ Saved: {h5_path}")     


def save_results_h5(index,x,y,path):

    h5_path = os.path.join(path, f"pressure_x{index}.h5")
    vtk_filename = f"parsed_vtk/pflotran_steady_pressure-001.vtk"
    with h5py.File(h5_path, "w") as h5f:
        h5f.create_dataset("grid points", data=list(zip(x,y)))   
        mesh = pv.read(vtk_filename)
        pressure = mesh.point_data['Liquid_Pressure']
        h5f.create_dataset('pressure', data = pressure)
        perm = mesh.point_data['Permeability']
        h5f.create_dataset('permeability', data = perm) 
    print(f" All results written to {h5_path}")