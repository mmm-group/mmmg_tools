import numpy as np
from mmmg_tools.objects import Structure
from pymatgen.core.structure import Molecule

# Read POSCAR
poscar = Structure.from_file('POSCAR')

# Create water molecule
x = 0.92 * np.sin(np.deg2rad(52.5))
z = 0.92 * np.cos(np.deg2rad(52.5))

molecule = Molecule(
        species = ['H','H','O'],
        coords = [[-x,0,z],[x,0,z],[0,0,0]],
        )

molecule.rotate_sites(theta=np.deg2rad(90 - poscar.structure.lattice.gamma), axis=[0,0,1])

# Select adsorption sites 
site = [
        poscar.structure._sites[128].coords,
        poscar.structure._sites[7].coords,
        ]

# Set adsorption c-vector 
angle = poscar.structure._sites[128].coords - poscar.structure._sites[259].coords

# Add the adsorbate
poscar.add_adsorbate(molecule, site[0], height=2.1, mirror_site=site[1], rotation=angle)

# Write POSCAR
poscar.write_file('add-POSCAR')
