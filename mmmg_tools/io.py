import numpy as np
from pymatgen.core import Lattice
from pymatgen.electronic_structure.bandstructure import BandStructure, Kpoint
from pymatgen.electronic_structure.core import Spin


def read_wavecar(filename):
    """
    Information is extracted from the given WAVECAR

    Args:
        filename (str): path to input file
    """

    # c = 0.26246582250210965422
    # 2m/hbar^2 in agreement with VASP
    _C = 0.262465832
    with open(filename, 'rb') as f:
        # read the header information
        recl, spin, rtag = np.fromfile(
            f, dtype=np.float64, count=3).astype(np.int)
        recl8 = int(recl / 8)

        # check to make sure we have precision correct
        if rtag != 45200 and rtag != 45210:
            raise ValueError('invalid rtag of {}'.format(rtag))

        # padding
        np.fromfile(f, dtype=np.float64, count=(recl8 - 3))

        # extract kpoint, bands, energy, and lattice information
        nk, nb, encut = np.fromfile(
            f, dtype=np.float64, count=3).astype(np.int)
        a = Lattice(np.fromfile(f, dtype=np.float64, count=9).reshape((3, 3)))
        efermi = np.fromfile(f, dtype=np.float64, count=1)[0]

        # calculate reciprocal lattice
        b = a.reciprocal_lattice._matrix

        # calculate maximum number of b vectors in each direction
        bmag = np.linalg.norm(b, axis=1)

        # calculate maximum integers in each direction for G
        phi12 = np.arccos(np.dot(b[0, :], b[1, :]) / (bmag[0] * bmag[1]))
        sphi123 = (np.dot(b[2, :], np.cross(b[0, :], b[1, :])) /
                   (bmag[2] * np.linalg.norm(np.cross(b[0, :], b[1, :]))))
        nbmaxA = np.sqrt(encut * _C) / bmag
        nbmaxA[0] /= np.abs(np.sin(phi12))
        nbmaxA[1] /= np.abs(np.sin(phi12))
        nbmaxA[2] /= np.abs(sphi123)
        nbmaxA += 1

        phi13 = np.arccos(np.dot(b[0, :], b[2, :]) / (bmag[0] * bmag[2]))
        sphi123 = (np.dot(b[1, :], np.cross(b[0, :], b[2, :])) /
                   (bmag[1] * np.linalg.norm(np.cross(b[0, :], b[2, :]))))
        nbmaxB = np.sqrt(encut * _C) / bmag
        nbmaxB[0] /= np.abs(np.sin(phi13))
        nbmaxB[1] /= np.abs(sphi123)
        nbmaxB[2] /= np.abs(np.sin(phi13))
        nbmaxB += 1

        phi23 = np.arccos(np.dot(b[1, :], b[2, :]) / (bmag[1] * bmag[2]))
        sphi123 = (np.dot(b[0, :], np.cross(b[1, :], b[2, :])) /
                   (bmag[0] * np.linalg.norm(np.cross(b[1, :], b[2, :]))))
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
            coeffs = [[[None for i in range(nb)]
                       for j in range(nk)] for _ in range(spin)]
            band_energy[Spin.down] = [
                [None for i in range(nk)] for j in range(nb)]
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
                enocc = np.fromfile(f, dtype=np.float64,
                                    count=3 * nb).reshape((nb, 3))
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
                    raise ValueError(
                        'failed to generate the correct number of G points')

                # extract coefficients
                for inb in range(nb):
                    if rtag == 45200:
                        data = np.fromfile(f, dtype=np.complex64, count=nplane)
                        buf = np.fromfile(f, dtype=np.float64,
                                          count=recl8 - nplane)
                    elif rtag == 45210:
                        # this should handle double precision coefficients
                        # but I don't have a WAVECAR to test it with
                        data = np.fromfile(
                            f, dtype=np.complex128, count=nplane)
                        np.fromfile(f, dtype=np.float64,
                                    count=recl8 - 2 * nplane)

                    if Gflag:
                        """
                        This is confusing but bear with me, occupations for most bands
                        are twice as large so need factor of sqrt(1/2). This isn't true
                        for the first band it seems this has some weird factor ~ the one
                        used.
                        """
                        data = data * \
                            np.sqrt(0.5) if inb != 0 else data * \
                            np.exp(-1 / (2 * np.pi))
                        data = np.concatenate((data, np.conj(data[1:])))
                        for iGp in range(len(Gpoints[ink])):
                            if iGp < len(Gptemp):
                                Gpoints[ink][iGp] = Gptemp[iGp].copy()
                            else:
                                Gpoints[ink][iGp] = -1 * \
                                    Gptemp[(iGp + 1) % len(Gptemp)].copy()

                    if spin == 2:
                        coeffs[ispin][ink][inb] = data
                    else:
                        coeffs[ink][inb] = data

    bands = BandStructure(kpoints, band_energy, a.reciprocal_lattice, efermi)

    return encut, bands, Gpoints, coeffs, _nbmax
