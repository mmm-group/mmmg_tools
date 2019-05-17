import numpy as np
from pymatgen.core import Lattice, ArrayWithUnit, PeriodicSite, Structure, FloatWithUnit
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
