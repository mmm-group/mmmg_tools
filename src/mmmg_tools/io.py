import numpy as np
from pymatgen.core import (
    Lattice, 
    ArrayWithUnit, 
    PeriodicSite, 
    Structure, 
    FloatWithUnit,
)
from pymatgen.electronic_structure.bandstructure import (
    BandStructure,
    Kpoint,
)
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import Poscar

def read_cube(filename):
    """
    Converts cube files to CHGCAR files.

    args:
        filename (str): Name of file.
    """
    comments, sites = [],[]
    with open(filename,'r') as f:
        nval=1
        while True:
            #Read until first proper line is printed
            l = f.readline().split()
            try:
                natoms=int(l[0])
                x0,y0,z0=[float(v) for v in l[1:4]]
                if len(l) > 4:
                    nval=float(l[5])
                break
            except ValueError:
                comments.append(' '.join(l))
        #Create Lattice object
        a = ArrayWithUnit([float(v) for v in f.readline().split()],'bohr')
        b = ArrayWithUnit([float(v) for v in f.readline().split()],'bohr')
        c = ArrayWithUnit([float(v) for v in f.readline().split()],'bohr')
        n1,n2,n3 = [int(v) for v in [a[0],b[0],c[0]]]
        npts = [n1,n2,n3]
        lattice = Lattice(
                [n1*a.to('ang')[1:],n2*b.to('ang')[1:],n3*c.to('ang')[1:]]
                )
        #Create Structure object
        for i in range(np.absolute(natoms)):
            l = f.readline().split()
            element = int(l[0])
            coord = ArrayWithUnit(
                    [float(v)+(0.5/npts[i]) for i,v in enumerate(l[-3:])],
                    'bohr')
            sites.append(PeriodicSite(
                (l[0]),coord.to('ang'),
                lattice,
                coords_are_cartesian=True))
        poscar = Poscar(Structure.from_sites(sites))
        #Check for DSET and collect m values
        if natoms < 0:
            l = f.readline().split()
            dset_ids = []
            while True:
                for m in l:
                    dset_id.append(int(m))
                if len(dset_id)+1 >= dset_id[0]:
                    break
                else:
                    l = f.readline().split()
            nval = dset_id[0]
        #Create CHG_data
        data_count,dc,nc = 0,0,0
        dataset = np.zeros([nval,n1,n2,n3])
        while dc < nval*n1*n2*n3:
            l = f.readline().split()
            for d in l:
                z = data_count % n3
                y = int(np.floor(data_count / n3)) % n2
                x = int(np.floor(data_count / n3 / n2))
                dataset[nc,x,y,z] = float(d)
                dc += 1
                nc = (nc + 1) % nval
                if nc == 0:
                    data_count += 1
        #Sum all orbitals if orbials decomposed
        if natoms < 0:
            for ds in dataset[1:]:
                dataset[0] = np.add(dataset[0],ds)
        data = {'total': dataset[0]}
        return (poscar, data, {})

def read_wavecar(filename='WAVECAR'):
    """
    Information is extracted from the given WAVECAR

    Args:
        filename (str): input file (default: WAVECAR)
        verbose (bool): determines whether processing information is shown
        precision (str): determines how fine the fft mesh is (normal or
                         accurate), only the first letter matters
    """

    # c = 0.26246582250210965422
    # 2m/hbar^2 in agreement with VASP
    _C = 0.262465832
    with open(filename, 'rb') as f:
        # read the header information
        recl, spin, rtag = np.fromfile(f, dtype=np.float64, count=3).astype(np.int)
        recl8 = int(recl / 8)

        # check to make sure we have precision correct
        if rtag != 45200 and rtag != 45210:
            raise ValueError('invalid rtag of {}'.format(rtag))

        # padding
        np.fromfile(f, dtype=np.float64, count=(recl8 - 3))

        # extract kpoint, bands, energy, and lattice information
        nk, nb, encut = np.fromfile(f, dtype=np.float64, count=3).astype(np.int)
        a = Lattice(np.fromfile(f, dtype=np.float64, count=9).reshape((3, 3)))
        efermi = np.fromfile(f, dtype=np.float64, count=1)[0]

        # calculate reciprocal lattice
        b = a.reciprocal_lattice._matrix

        # calculate maximum number of b vectors in each direction
        bmag = np.linalg.norm(b, axis=1)

        # calculate maximum integers in each direction for G
        phi12 = np.arccos(np.dot(b[0,:], b[1,:]) / (bmag[0] * bmag[1]))
        sphi123 = (np.dot(b[2, :], np.cross(b[0,:], b[1,:])) / 
            (bmag[2] * np.linalg.norm(np.cross(b[0,:], b[1,:]))))
        nbmaxA = np.sqrt(encut * _C) / bmag
        nbmaxA[0] /= np.abs(np.sin(phi12))
        nbmaxA[1] /= np.abs(np.sin(phi12))
        nbmaxA[2] /= np.abs(sphi123)
        nbmaxA += 1

        phi13 = np.arccos(np.dot(b[0,:], b[2,:]) / (bmag[0] * bmag[2]))
        sphi123 = (np.dot(b[1,:], np.cross(b[0,:], b[2,:])) /
            (bmag[1] * np.linalg.norm(np.cross(b[0,:], b[2,:]))))
        nbmaxB = np.sqrt(encut * _C) / bmag
        nbmaxB[0] /= np.abs(np.sin(phi13))
        nbmaxB[1] /= np.abs(sphi123)
        nbmaxB[2] /= np.abs(np.sin(phi13))
        nbmaxB += 1

        phi23 = np.arccos(np.dot(b[1, :], b[2, :]) / (bmag[1] * bmag[2]))
        sphi123 = (np.dot(b[0,:], np.cross(b[1,:], b[2,:])) / 
            (bmag[0] * np.linalg.norm(np.cross(b[1,:], b[2,:]))))
        nbmaxC = np.sqrt(encut * _C) / bmag
        nbmaxC[0] /= np.abs(sphi123)
        nbmaxC[1] /= np.abs(np.sin(phi23))
        nbmaxC[2] /= np.abs(np.sin(phi23))
        nbmaxC += 1

        _nbmax = np.max([nbmaxA, nbmaxB, nbmaxC], axis=0).astype(np.int)

        # padding
        np.fromfile(f, dtype=np.float64, count=recl8 - 13)

        # reading records
        # np.set_printoptions(precision=7, suppress=True)
        Gpoints = [None for _ in range(nk)]
        kpoints = []
        band_energy = {Spin.up: [[None for i in range(nk)] for j in range(nb)]}
        if spin == 2:
            coeffs = [[[None for i in range(nb)] for j in range(nk)] for _ in range(spin)]
            band_energy[Spin.down] = [[None for i in range(nk)] for j in range(nb)]
        else:
            coeffs = [[None for i in range(nb)] for j in range(nk)]

        for ispin in range(spin):
            for ink in range(nk):
                # information for this kpoint
                nplane = int(np.fromfile(f, dtype=np.float64, count=1)[0])
                kpoint = np.fromfile(f, dtype=np.float64, count=3)

                if ispin == 0:
                    kpoints.append(kpoint)
                else:
                    assert np.allclose(kpoints[ink], kpoint)

                # energy and occupation information
                enocc = np.fromfile(f, dtype=np.float64, count=3 * nb).reshape((nb, 3))
                if ispin == 1:
                    for iband, eig in enumerate(enocc):
                        band_energy[Spin.down][iband][ink] = eig[0]
                else:
                    for iband, eig in enumerate(enocc):
                        band_energy[Spin.up][iband][ink] = eig[0]

                # padding
                np.fromfile(f, dtype=np.float64, count=(recl8 - 4 - 3 * nb))

                # generate G integers
                gpoints = []
                for i in range(2 * _nbmax[2] + 1):
                    i3 = i - 2 * _nbmax[2] - 1 if i > _nbmax[2] else i
                    for j in range(2 * _nbmax[1] + 1):
                        j2 = j - 2 * _nbmax[1] - 1 if j > _nbmax[1] else j
                        for k in range(2 * _nbmax[0] + 1):
                            k1 = k - 2 * _nbmax[0] - 1 if k > _nbmax[0] else k
                            G = np.array([k1, j2, i3])
                            v = kpoint + G
                            g = np.linalg.norm(np.dot(v, b))
                            E = g ** 2 / _C
                            if E < encut:
                                gpoints.append(G)
                Gpoints[ink] = np.array(gpoints, dtype=np.float64)
                Gflag = False
                if (len(Gpoints[ink]) + 1)/2 == nplane:
                    Gflag = True
                    gptemp = []
                    for i in range(2 * _nbmax[2] + 1):
                        i3 = i - 2 * _nbmax[2] - 1 if i > _nbmax[2] else i
                        for j in range(2 * _nbmax[1] + 1):
                            j2 = j - 2 * _nbmax[1] - 1 if j > _nbmax[1] else j
                            for k in range(_nbmax[0] + 1):
                                if k == 0 and j2 < 0:
                                    pass
                                elif k == 0 and j2 == 0 and i3 < 0:
                                    pass
                                else:
                                    G = np.array([k, j2, i3])
                                    v = kpoint + G
                                    g = np.linalg.norm(np.dot(v, b))
                                    E = g ** 2 / _C
                                    if E < encut:
                                        gptemp.append(G)
                    Gptemp = np.array(gptemp, dtype=np.float64)
                elif len(Gpoints[ink]) != nplane:
                    print(len(Gpoint), nplane)
                    raise ValueError('failed to generate the correct number of G points')

                # extract coefficients
                for inb in range(nb):
                    if rtag == 45200:
                        data = np.fromfile(f, dtype=np.complex64, count=nplane)
                        buf = np.fromfile(f, dtype=np.float64, count=recl8 - nplane)
                    elif rtag == 45210:
                        # this should handle double precision coefficients
                        # but I don't have a WAVECAR to test it with
                        data = np.fromfile(f, dtype=np.complex128, count=nplane)
                        np.fromfile(f, dtype=np.float64, count=recl8 - 2 * nplane)

                    if Gflag:
                        data = data * np.sqrt(0.5) if inb != 0 else data * np.exp(-1 / (2 * np.pi))
                        data = np.concatenate((data, np.conj(data[1:])))
                        for iGp in range(len(Gpoints[ink])):
                            if iGp < len(Gptemp):
                                Gpoints[ink][iGp] = Gptemp[iGp].copy()
                            else:
                                Gpoints[ink][iGp] = -1 * Gptemp[(iGp + 1) % len(Gptemp)].copy()

                    if spin == 2:
                        coeffs[ispin][ink][inb] = data
                    else:
                        coeffs[ink][inb] = data

    bands = BandStructure(kpoints, band_energy, a.reciprocal_lattice, efermi) 

    return encut, bands, Gpoints, coeffs, _nbmax
