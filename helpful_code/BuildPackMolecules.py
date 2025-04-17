#create boxes of packed moelcuels for simulation runs. Note that this does not equilibrate the structures, just creates the top and gro files based on experimental densityies with some packing room 

from openff.toolkit import Molecule, Topology, ForceField
from openff.interchange import Interchange
from openff.units import unit, Quantity 
import numpy as np

import openff.nagl
from openff.nagl import GNNModel 
from openff.nagl_models import list_available_nagl_models
from openff.interchange.components._packmol import pack_box


# Make molecules
methylaetate = Molecule.from_smiles("CC(=O)OC", allow_undefined_stereo=True)


alkanes = {
    'methylacetate': methylaetate
}

# Experimental densities in g/nm^3 for approximate packing
densities = {
    'methylacetate': 9.3e-22 / 74.08, #g/nm^3 divided by g/mol to get mol/nm^3
    # 'hexane': 6.59e-22 / 86.17, 
    # 'heptane': 6.8e-22 / 100.21, 
    # 'octane': 7.03e-22 / 114.23, 
    # 'decane': 7.3e-22 / 152.29, 
    # 'pentadecane': 7.69e-22 / 212.42, 
}

size = 1000  # Number of molecules

# Dictionary to store single box size values for liquid simulations
box_sizes = {}

# Calculate box sizes for each molecule
for molecule, density in densities.items():
    box_size = ((size / (density * 6.02e23)) * 2) ** (1/3)  # Add packing room. 
    box_sizes[molecule] = round(box_size, 2)  # Store as single float

# for molecule, box_size in box_sizes.items():
#     print(f"Box size for {molecule}: {box_size} nm")

for name, alkane in alkanes.items():
    alkane.generate_conformers(n_conformers=1)
    # print(f"{name} conformers: {len(alkane.conformers)}")

    if not alkane.conformers:
        print(f"Error: No conformers generated for {name}")
        continue

    alkane.name = "MEA"
    for i, atom in enumerate(alkane.atoms, 3):
        atom.metadata["residue_name"] = 'MEA'
    alkane.generate_unique_atom_names()

    # Get the gas phase. Lets use a 10 nm box 
    gas_box = unit.Quantity(100 * np.eye(3), unit.angstrom)
    gas_top = Topology.from_molecules([alkane])
    ff = ForceField("openff-2.2.0.offxml")
    gas_interchange = Interchange.from_smirnoff(
        force_field=ff,
        topology=gas_top,
        box=gas_box,
    )
    gas_interchange.to_gromacs(f"{name}_gas")

    # Get the liquid gromacs inputs
    cubic_box_size = box_sizes[name]  # Directly access float
    
    cubic_box = unit.Quantity(10 * cubic_box_size * np.eye(3), unit.angstrom)  # Convert nm to Å


    try:
        packed_topol = pack_box(
            molecules=[alkane], number_of_copies=[size], solute=None,
            tolerance=2 * unit.angstrom, box_vectors=cubic_box
        )
        packed_interchange = Interchange.from_smirnoff(
            force_field=ff,
            topology=packed_topol,
            box=cubic_box
        )
        packed_interchange.to_gromacs(f"{name}_liquid")
    except Exception as e:
        print(f"Error packing liquid box for {name}: {e}")
