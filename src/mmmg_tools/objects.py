import numpy as np
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
            unit = poscar.structure.lattice.volume / FloatWithUnit(1, 'borh^3').to('ang^3')
            for k in data.keys():
                data[k] *= unit
        elif ftype.lower() == 'chgcar':
            poscar, data, data_aug = vd.parse_file(filename)
        return cls(poscar, data, data_aug=data_aug)

    @classmethod
    def from_wavefunction(cls, filename, ftype='WAVECAR', structure=None, bands=None, erange=None):
        """
        Create Charge object from wavefunctions. 
        
        See pymatgen.io.vasp.outputs.Wavecar.get_parchg() for reasoning as to why the charge is 
        incorrect.
        
        args:
            filename (str): path to wavefunction file.
            ftype (str): see Wavefunction class doc for supported file-types. Defualt = WAVECAR.
            structure (Structure): associated structure object if required.
            band (List): list of bands to include by band index. Defualt = All <= Ef.
            erange (list): list containing upper and lower bounds of included band energies
                           relative to Ef. Default = All <= Ef.
        """
        if bands is not None and erange is not None:
            print("Both bands and erange are set, only using bands.")
        wav = Wavefunction.from_file(filename, ftype=ftype, structure=structure)
        data = {'total': np.zeros(tuple(wav.ng*2))}
        if bands is not None:
            data['diff'] = np.zeros(tuple(wav.ng*2))
            chg = cc(wav.structure, data)
            for ik in range(wav.nk):
                for ib in bands:
                    chg += wav.get_parchg(ik, ib)
            return cls(wav.structure, chg.data.copy())
        elif erange is not None:
            emin = erange[0] + wav.efermi
            emax = erange[1] + wav.efermi
        else:
            emin = -float('inf')
            emax = wav.efermi
        chg_up = cc(wav.structure, data.copy())
        if wav.spin == 2:
            chg_down = cc(wav.structure, data.copy())
            for ik in range(wav.nk):
                for ib in range(wav.nb):
                    if emin <= wav.band_energy[Spin.up][ib][ik] <= emax:
                        chg_up += wav.get_parchg(ik, ib, spin=0)
                    elif wav.band_energy[Spin.up][ib][ik] > emax:
                        break
                    if emin <= wav.band_energy[Spin.down][ib][ik] <= emax:
                            chg_down += wav.get_parchg(ik, ib, spin=1)
                    elif wav.band_energy[Spin.down][ib][ik] > emax:
                        break
            data['total'] = chg_up.data['total'] + chg_down.data['total']
            data['diff'] = chg_up.data['total'] - chg_down.data['total']
        else:
            for ik in range(wav.nk):
                for ib in range(wav.nb):
                    if emin <= wav.band_energy[Spin.up][ib][ik] <= emax:
                        chg_up += wav.get_parchg(ik, ib, spin=0)
            data = chg_up.data
        return cls(wav.structure, data)

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
        Read in data from a potential file (LOCPOT).

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

    def get_parchg(self, kpoint, band, spin=None, phase=False, scale=2):
        """
        Wrapper for pymatgen function of same name changing the returned object.

        Args:
            kpoint (int):   the index of the kpoint for the wavefunction
            band (int):     the index of the band for the wavefunction
            spin (int):     optional argument to specify the spin. If the
                            Wavecar has ISPIN = 2, spin == None generates a
                            chgcar with total spin and magnetization, and
                            spin == {0, 1} specifies just the spin up or
                            down component.
            phase (bool):   flag to determine if the charge density is
                            multiplied by the sign of the wavefunction.
                            Only valid for real wavefunctions.
            scale (int):    scaling for the FFT grid. The default value of 2 is
                            at least as fine as the VASP default.

        Returns:
            Charge object
        """

        _chg = super().get_parchg(self.structure, kpoint, band, spin=spin, phase=phase, scale=scale)
        return Charge(self.structure, _chg.data.copy())
