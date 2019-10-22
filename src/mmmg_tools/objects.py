import numpy as np
from scipy.signal import resample
import pandas as pd
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import cpu_count
from tqdm.auto import tqdm
from pymatgen.io.vasp.outputs import (
    Chgcar as cc, 
    Locpot as lp,
    VolumetricData as vd,
    Wavecar as wc,
)
from pymatgen.io.vasp.inputs import Poscar as pc
from pymatgen.core import (
    Site, 
    FloatWithUnit, 
    Structure as st, 
    Lattice, 
    PeriodicSite,
)
from pymatgen.electronic_structure.core import Spin
from .io import (
    read_cube,
    read_wavecar,
)
from .interfaces import Bader
from .pybader import bader_wrap as bw
from functools import partial

class Charge(cc):
    """
    Charge file object.

    Supported file-types:
        - CHGCAR
        - cube
    """

    def __init__(self, structure, data, data_aug=None, voxel_origin=(0.,)*3):
        self._voxel_origin = voxel_origin
        self._voxel_lattice = Lattice(
            np.divide(structure.structure.lattice._matrix, data['total'].shape)
        )
        cc.__init__(self, Structure.from_structure(structure.structure), data, data_aug)

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
            unit = poscar.structure.lattice.volume / FloatWithUnit(1, 'bohr^3').to('ang^3')
            data = {k: v * unit for k, v in data.items()}
            vox = (.5,) * 3
        elif ftype.lower() == 'chgcar':
            poscar, data, data_aug = vd.parse_file(filename)
            vox = (0.,) * 3
        return cls(poscar, data, data_aug=data_aug, voxel_origin=vox)

    @classmethod
    def from_wav(cls, wavefunction, kweight, blist=None, brange=None, klist=None, spin=None):
        """
        Create Charge object from wavefunctions. 
        
        See pymatgen.io.vasp.outputs.Wavecar.get_parchg() for reasoning as to
        why the charge is incorrect.
        
        args:
            wavefunction (object): Wavefunction object.
            kweight (List): list of kpoint weights.
            blist (List):   list of bands to include by band index. 
                            Defualt = All <= Ef.
            brange (List):  list containing upper and lower bounds of included 
                            band energies relative to Ef. Default = All <= Ef.
            klist (List):   list of kpoints to include. Default = All.
            spin (int):     0 or 1 to choose spin up or down. Default Both.

        """
        sdict, kdict = {}, {}
        if klist is None:
            klist = range(wavefunction.nk)
        spin = list([0, 1] if spin is None else [spin])
        if blist is not None:
            for s in spin:
                for ik in klist:
                    kdict[ik] = blist
                sdict[s] = kdict
        elif brange is not None:
            emin, emax = tuple(np.add(brange, wavefunction.efermi))
            for s in spin:
                sp = Spin.up if s == 0 else Spin.down
                for ik in klist:
                    kdict[ik] = [
                        ib for ib in range(wavefunction.nb) 
                        if emin <= wavefunction.band_energy[sp][ib][ik] <= emax
                    ]
                sdict[s] = kdict
        else:
            emin = -float('inf')
            emax = 0
            for s in spin:
                sp = Spin.up if s == 0 else Spin.down
                for ik in klist:
                    kdict[ik] = [
                        ib for ib in range(wavefunction.nb) 
                        if emin <= wavefunction.band_energy[sp][ib][ik] <= emax
                    ]
                sdict[s] = kdict
        data = wavefunction.get_density(sdict, kweight)
        return cls(wavefunction.structure, data, voxel_origin=(0.)*3)

    @property
    def voxel_origin(self):
        return self._voxel_origin

    @property
    def voxel_lattice(self):
        return self._voxel_lattice

    def bader_calc(self, rho=None, ref=None, args=None):
        """
        Run Bader charge analysis on a charge density. 

        args:
            rho (ndarray): Charge density. Default = self.data['total'].
            ref (ndarray): Reference density. Default = self.data['total'].
            args (str): Command line arguments. Default = None.
        """
        lat = self.structure.lattice.matrix
        ion_shift = np.divide(self.voxel_origin, self.dim)
        ions = np.add(self.structure.frac_coords, ion_shift)
        self.bader = Bader(self.structure, self.dim)
        if rho is None:
            rho = self.data['total']
        if ref is None:
            ref = rho
        if args is None:
            args = ''
        bw.bader_run(lat, rho, ref, ions, args, self.bader, self.bader.masks)

    def to_structure(self):
        """
        Produces Structure object.
        """
        return Structure.from_structure(self.structure)

class Structure(pc):
    """
    Structure object.

    Supported file-types:
        - POSCAR
    """

    @classmethod
    def from_structure(cls, structure):
        """
        Create Structure object from pymatgen.core.Structure.

        args:
            structure (object): Pymatgen structure object.
        """
        return cls(structure)

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
            ad_pb_atoms (bool): whether to append the boundary atoms. Default = False
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

    def sort(self, reverse=False):
        """
        Organise sites list based on site species.

        args:
            reverse (bool): sorts reverse alphabetically. Default = False
        """
        self.structure._sites.sort(key=lambda site: site.specie, reverse=reverse) 

    def generalised_coordination_number(self, order:int=1, bond_cut:float=2.5):
        """
        Get the genealised cooordition number as presented in:
        https://onlinelibrary.wiley.com/doi/10.1002/anie.201402958
        
        Requires c_max as a site property: 
        self.structure.add_site_property('c_max', [6 for s in self.structure.sites])

        args:
            order: the order of the correction to the coordination number.
            bond_cut: the distance under which to consider an atom coordinated.
        """
        coord = list()
        neighbour = list()
        for site in self.structure.sites:
            neighbour.append(self.structure.get_neighbors(site, bond_cut, include_index=True))
            coord.append(len(neighbour[-1]))
# py3.8 while i:=0 < order
        i = 0
        while i < order:
            new_coord = list()
            for j, site in enumerate(self.structure.sites):
                c = list()
                for site in neighbour[j]:
                    c.append(coord[site[2]] / site.c_max)
                new_coord.append(np.sum(c))
            coord = new_coord.copy()
            i += 1
        self.structure.add_site_property('GCN', coord)

    def write_xyz(self, filename:str, site_properties:list=[]):
        """
        Write Structure as an xyz.

        args:
            filename: where to save file to.
            site_properties: list of properties to append to file as extra columns.
        """
        string = f"{np.sum(self.natoms)}\n{'  '.join(set(self.site_symbols))}"
        for site in self.structure.sites:
            string += f"\n{site.specie}  {site.x}  {site.y}  {site.z}"
            for prop in site_properties:
                string += f"  {site.__dict__['properties'][prop]}"
        with open(filename, 'wt') as f:
            f.write(string)

class Potential(lp):
    """
    Potential file object.

    Supported file-types:
        - LOCPOT
        - cube
    """

    def __init__(self, structure, data, voxel_origin=(0.,)*3):
        self._voxel_origin = voxel_origin
        self._voxel_lattice = Lattice(
            np.divide(structure.structure.lattice._matrix, data['total'].shape)
        )
        lp.__init__(self, Structure.from_structure(structure.structure), data)

    @classmethod
    def from_file(cls, filename, ftype='LOCPOT'):
        """
        Read in data from a potential file (LOCPOT).

        args:
            filename (str): path to file.
            ftype (str): See class doc for supported file-types. Default = LOCPOT.
        """
        if ftype.lower() == 'cube' or filename.split('.')[-1] == 'cube':
            poscar, data, data_aug = read_cube(filename)
            unit = FloatWithUnit(1, 'Ha').to('eV') 
            data = {k: v * unit for k, v in data.items()}
            vox = (.5,) * 3
        elif ftype.lower() == 'locpot':
            poscar, data, data_aug = vd.parse_file(filename)
            vox = (0.,) * 3
        return cls(poscar, data, voxel_origin=vox)

    @property
    def voxel_origin(self):
        return self._voxel_origin

    @property
    def voxel_lattice(self):
        return self._voxel_lattice

    def to_structure(self):
        """
        Produces Structure object.
        """
        return Structure.from_structure(self.structure)

    def onsite_electrostatic(self, Z, r=1.2, sites=None, upsample=1):
        """
        Calculates the on-site potential for atoms within the structure.

        args:
            Z (float/dict): normalisation for s-orbital, can be supplied as dict for species.
            r (float/dict): radius ofsphere around atom, can a dict for species.
            sites (list): list of int or sites to get onsite potential of.
            upsample (int): scale to upsample the data by using scipy.resample.  
        """
        if sites is None:
            sites = self.structure.sites
        else:
            sites = [self.structure._sites[s] if isinstance(s, int) else s for s in sites]
        if isinstance(Z, float):
            _set = set([s.specie.symbol for s in sites])
            Z = {k: Z for k in _set}
        if isinstance(r, float):
            _set = set([s.specie.symbol for s in sites])
            r = {k: r for k in _set}
        data = self.data['total'].copy()
        for i, dim in enumerate(np.multiply(self.dim, upsample)):
            data = resample(data, dim, axis=i)
        onsite = list()
        specie = list()
        lat = self.voxel_lattice.scale(self.voxel_lattice.volume / upsample**3)
        for site in tqdm(sites):
            specie.append(site.specie.symbol)
            z = Z[specie[-1]]
            _c = lat.get_fractional_coords(site.to_unit_cell().coords)
            coords = np.array(_c) - np.array(_c, dtype=int)
            _,dist,_, image = lat.get_points_in_sphere(
                [self.voxel_origin],
                lat.get_cartesian_coords(coords),
                r[specie[-1]],
                zip_results=False,
            )
            idx = np.mod(np.add(np.array(_c, dtype=int), image), np.multiply(self.dim, upsample))
            poten = data[idx[:,0], idx[:,1], idx[:,2]]
            onsite.append(np.multiply(np.exp(-2 * z * dist), poten).sum())
        columns = ['Species', 'Electrostatic Potential']
        index = [self.structure._sites.index(site) for site in sites]
        onsite = np.multiply(onsite, z**3 * lat.volume / np.pi)
        data = zip(specie, onsite)
        return pd.DataFrame(data, index, columns)

class Wavefunction(wc):
    """
    Wave-function object.

    Supported file-types:
        - WAVECAR (requires Structure)
    """

    def __init__(self, encut, bands, Gpoints, coeffs, _nbmax, structure):
        self._C = 0.262465832
        self._bands = bands
        self._lattice = bands.lattice_rec
        self._structure = structure
        self.spin = 2 if bands.is_spin_polarized else 1
        self.nk = len(bands.kpoints)
        self.nb = bands.nb_bands
        self.encut = encut 
        self.efermi = bands.efermi
        self.a = self._lattice.reciprocal_lattice._matrix
        self.b = self._lattice._matrix
        self.kpoints = bands.kpoints
        self.band_energy = bands.bands
        self.Gpoints = Gpoints
        self.coeffs = coeffs
        self._nbmax = _nbmax
        self.ng = 4 * _nbmax

    @property
    def structure(self):
        return self._structure

    @classmethod
    def from_file(cls, filename, ftype='WAVECAR', structure=None):
        """
        Read in wavefunctions from file. 
        
        Some file formats require a structure object to be passed also as the file
        doesn't contain this information.

        args:
            filename (str): path to file.
            ftype (str): See class doc for supported file-types. Default = WAVECAR.
            structure (Structure): a structure object if required by file-type.
        """
        if ftype.lower() == 'wavecar':
            if structure is None:
                raise ValueError('structure arg is required for WAVECAR file-types.')
            encut, bands, Gpoints, coeffs, _nbmax = read_wavecar(filename)
        return cls(encut, bands, Gpoints, coeffs, _nbmax, structure)

    def get_density(self, sdict, kweight, scale=2):
        """
        Produce density for Charge object.
        """
        data, den= {}, {}
        ng = tuple(self.ng * scale)
        Ng = np.prod(ng)
        for spin, kdict in sdict.items():
            if all(s == 1 for s in [spin, self.spin]):
                return data
            den[spin] = np.zeros(ng)
            origin = [l // 2 for l in ng]
            print(f"Spin {Spin(np.sign((-2*spin)+1))._name_}: ")
            for ik, nb in tqdm(kdict.items()):
                coeffs = np.array(
                    self.coeffs[spin][ik] if self.spin == 2 else self.coeffs[ik]
                )
                gpoint = np.add(origin, self.Gpoints[ik].astype(np.int))
                def wavgen(ib):
                    mesh = np.zeros(ng, dtype=np.complex)
                    for g, c in zip(gpoint, coeffs[ib]):
                        mesh[tuple(g)] = c
                    wf = np.fft.ifftn(mesh, norm='ortho')
                    wf /= np.sqrt(np.vdot(wf, wf))
                    return (np.conj(wf) * wf).real
                pool = Pool(cpu_count() - 1)
                rho = pool.map(wavgen, nb)
                pool.close()
                den[spin] +=  np.multiply(np.sum(rho, axis=0), kweight[ik])
        if len(sdict) == 2:
            data['total'] = den[0] + den[1]
            data['diff'] = den[0] - den[1]
        else:
            data['total'] = den[spin]
        return {k: np.multiply(v, Ng) for k, v in data.items()}

    def real_space(self, coeffs, gpoints):
        ng = tuple(self.ng * 2)
        mesh = np.zeros(ng, dtype=np.complex)
        for g, c in zip(map(tuple, gpoints), coeffs):
            mesh[g] = c
        wfr = np.fft.ifftn(mesh, norm='ortho')
        wfr /= np.sqrt(np.vdot(wfr, wfr))
        return wfr
