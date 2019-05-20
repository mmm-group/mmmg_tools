import numpy as np
from pymatgen.io.vasp.outputs import Chgcar as cc, Locpot as lp
from pymatgen.io.vasp.outputs import VolumetricData as vd
from pymatgen.io.vasp.inputs import Poscar as pc
from pymatgen.core import Site, FloatWithUnit, Structure as st, Lattice, PeriodicSite
from .io import read_cube
from .interfaces import Bader
from .pybader import bader_wrap as bw

class Charge(cc):
    """
    Charge file object.

    Supported file-types:
        - CHGCAR
        - cube
    """

    @classmethod
    def from_file(cls, filename, ftype='CHGCAR'):
        """
        Read in data from chgcar or cube file.

        args:
            filename (str): path to file.
            ftype (str): See class doc for supported file-types. Default = CHGCAR.
        """
        if ftype.lower() == 'cube' or filename.split('.')[-1] == 'cube':
            poscar, data, data_aug = read_cube(filename)
            data = {k: v / 
                    (FloatWithUnit(1,'bohr^3').to('ang^3') / 
                        poscar.structure.lattice.volume) for k, v in data.items()}
        elif ftype.lower() == 'chgcar':
            poscar, data, data_aug = vd.parse_file(filename)
        return cls(poscar, data, data_aug=data_aug)

    def bader_calc(self, rho=None, ref=None, args=None):
        """
        Run Bader charge analysis on a charge density. 

        args:
            rho (ndarray): Charge density. Default = self.data['total'].
            ref (ndarray): Reference density. Default = self.data['total'].
            args (str): Command line arguments. Default = None.
        """
        lat = self.structure.lattice.matrix
        ions = self.structure.frac_coords
        self.bader = Bader(self.structure, self.dim)
        if rho is None:
            rho = self.data['total']
        if ref is None:
            ref = rho
        if args is None:
            args = ''
        bw.bader_run(lat, rho, ref, ions, args, self.bader, self.bader.masks)

class Structure(pc):
    """
    Structure object.

    Supported file-types:
        - POSCAR
    """

    @classmethod
    def from_file(cls, filename, ftype='POSCAR'):
        """
        Read in structure from file.

        args:
            filename (str): path to file.
            ftype (str): See class doc for supported file-types. Default = POSCAR.
        """
        if ftype.lower() == 'poscar':
            with open(filename,'r') as f:
                poscar = pc.from_string(f.read())
            data = poscar.as_dict()
        return cls(st.from_dict(data['structure']),
                comment = data['comment'],
                selective_dynamics = data['selective_dynamics'],
                true_names = data['true_names'],
                velocities = data.get("velocities", None),
                predictor_corrector = data.get("predictor_corrector", None),
                )

    def add_adsorbate(self, molecule, site, height=0, mirror_site=None, rotation=None):
        """
        Adosorb molecule on surface.

        args:
            molecule (Molecule): The adsorbate as a pymatgen object.
            site (3x1 array): Position to adsorb in cartesian coordinates.
            height (float): Addition z-height from site. Default = 0.
            mirror_site (3x1 array): Cartesian coords of site. Default = None.
            rotation (3x1 array): Vector to rotate the z-axis too. Default = None.
        """
        if rotation is not None:
            r = rotation/np.linalg.norm(rotation)
            theta, b = np.arccos(r[2]), np.cross([0,0,1],r)
            b = b/np.linalg.norm(b)
            c, s = np.cos(theta/2), np.sin(theta/2)
            q = [c,*np.multiply(s,b)]
            Q = np.array([
                [
                    q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 
                    2*(q[1]*q[2] - q[0]*q[3]), 
                    2*(q[1]*q[3] + q[0]*q[2]),
                ],
                [
                    2*(q[1]*q[2] + q[0]*q[3]), 
                    q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2, 
                    2*(q[2]*q[3] - q[0]*q[1]),
                ],
                [
                    2*(q[1]*q[3] - q[0]*q[2]),
                    2*(q[2]*q[3] + q[0]*q[1]),
                    q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2,
                ],
            ])
        else:
            Q = np.array([[1,0,0],[0,1,0],[0,0,1]])

        if mirror_site is None:
            for s in molecule:
                coords = np.dot(Q, np.add([0,0,height], s.coords).T)
                self.structure.append(
                    s.specie,
                    coords + site,
                    coords_are_cartesian=True,
                    properties=s.properties)
        else:
            for s in molecule:
                coords = np.dot(Q, np.add([0,0,height], s.coords).T)
                self.structure.append(
                    s.specie,
                    coords + site,
                    coords_are_cartesian=True,
                    properties=s.properties)
                self.structure.append(
                    s.specie,
                    mirror_site - coords,
                    coords_are_cartesian=True,
                    properties=s.properties)

    def set_vacuum(self, vacuum, add_pb_atoms=False):
        """
        Sets the vacuum in the z-axis of a Structure.

        args:
            vacuum (float): Amount of vacuum to add.
            ad_pb_atoms (bool): wether to append the boundary atoms. Default = False
        """

        if add_pb_atoms: 
            gen = (s for s in self.structure.sites 
                    if np.round(s.z) in (0, np.round(s.lattice.c)))
            for s in gen:
                coord = s.frac_coords
                coord[2] += 1 if s.c < .5 else -1
                self.structure.append(s.specie, coord, properties=s.properties)

        l = float('inf')
        h = .0

        for s in self.structure.sites:
            l = s.z if s.z < l else l
            h = s.z if s.z > h else h

        s = self.structure
        l = Lattice.from_lengths_and_angles(
                abc = [s.lattice.a, s.lattice.b, (h - l + vacuum)],
                ang = s.lattice.angles)
        new_sites = []
        for site in s.sites:
            new_site = PeriodicSite(site.species, site.coords, l,
                    properties=site.properties, coords_are_cartesian=True,
                    to_unit_cell=False)
            new_sites.append(new_site)
        self.structure._sites = new_sites
        self.structure._lattice = l

class Potential(lp):
    """
    Potential file object.

    Supported file-types:
        - LOCPOT
        - cube
    """

    @classmethod
    def from_file(cls, filename, ftype='LOCPOT'):
        """
        Read in data from locpot or cube file.

        args:
            filename (str): path to file.
            ftype (str): See class doc for supported file-types. Default = LOCPOT.
        """
        if ftype.lower() == 'cube' or filename.split('.')[-1] == 'cube':
            poscar, data, data_aug = read_cube(filename)
            data = {k: v / FloatWithUnit(1,'bohr^3').to('ang^3') for k, v in data.items()}
        elif ftype.lower() == 'locpot':
            poscar, data, data_aug = vd.parse_file(filename)
        return cls(poscar, data)
