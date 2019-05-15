import numpy as np
import pandas as pd
from .utils import sparray

class Bader():
    """
    Bader callback class.

    properties:
        mask: Mask class.
        bader_charge (list): Charge of each Bader volume.
        bader_ion (list): Associated ion for each Bader volume.
        bader_distance (list): Distance from each Bader maxima to ion.
        bader_pos (list): Position of each Bader maxima.
        ion_charge (list): Charge of each Bader ion.
        ion_volume (list): Volume of each Bader ion.
        ion_minsurfdist (list): Minimum distance to the Bader surface.
        vacuum_charge (float): Charge of the vacuum.
        vacuum_volume (float): Volume of the vacuum.
    """
    class Mask():
        """
        Mask callback class.

        properties:
            bader (list): List of the masks for each Bader volume.
        """
        def __init__(self,dim):
            """
            Initialise the class and bader property.

            args:
                dim (3x1 array): Dimensions to initialise sparse arrays.
            """
            self.dim = dim
            self.bader = []

        def __call__(self,vol,mask):
            """
            The callback function.
            """
            self.bader.append(sparray(self.dim,data=vol))

        def ion_mask(self,site,ion_map):
            """
            Return a sparse matrix representation of the ion mask.

            The ion mask is the array with which to multiply the charge
            density by. This then creates the charge density of ion 
            which can be filled and then written to file.
            """
            imask = sparray(self.dim)
            gen = (vol for vol,ion in enumerate(ion_map) if ion == site)
            for vol in gen:
                imask +=  self.bader[vol]
            return imask

    def __init__(self,structure,mask_dim):
        """
        Initialise the class and call the Mask class.

        args:
            structure (Structure): A pymatgen Structure object.
            mask_dim (3x1 array): Dimensions to intialise sparse array.
        """
        self.structure = structure
        self.masks = self.Mask(mask_dim)

    def __call__(self,nb,ni,bc,bi,bdist,bpos,ic,iv,imsd,vc,vv):
        """
        The callback function.
        """
        self.bader_charge = bc.copy()
        self.bader_ion = bi.copy()
        self.bader_distance = bdist.copy()
        self.bader_pos = bpos.copy()
        self.ion_charge = ic.copy()
        self.ion_volume = iv.copy()
        self.ion_minsurfdist = imsd.copy()
        self.vacuum_charge = vc
        self.vacuum_volume = vv
        data = {
                'Element' : self.structure.species,
                'X' : self.structure.frac_coords[:,0],
                'Y' : self.structure.frac_coords[:,1],
                'Z' : self.structure.frac_coords[:,2],
                'Charge' : self.ion_charge,
                'Min Surface Dist' : self.ion_minsurfdist,
                'Volume' : self.ion_volume,
                }
        self.ion_df = pd.DataFrame(data=data)
        data = {
                'X' : self.bader_pos[:,0],
                'Y' : self.bader_pos[:,1],
                'Z' : self.bader_pos[:,2],
                'Charge' : self.bader_charge,
                'Atom' : self.bader_ion,
                'Distance' : self.bader_distance,
                }
        self.bader_df = pd.DataFrame(data=data)
        del(self.structure,data)

    @property
    def nelect(self):
        """
        Return the total number of electrons.
        """
        return self.ion_df['Charge'].sum()

    @property
    def ACF(self):
        """
        Display the equivalent of the ACF.dat file.
        """
        idx = pd.Series(range(1,len(self.ion_df.index)+1))
        _df = self.ion_df.iloc[:,1:].rename_axis(columns='#')
        _df = _df.set_index(idx)
        txt = _df.to_string(justify='center').split('\n')
        for i,t in enumerate(txt):
            txt[i] = ' ' + t
        txt.insert(1,'-'*(len(txt[1])+1))
        txt.append('-'*len(txt[1]))
        txt.append(f"   Vacuum Charge: {self.vacuum_charge:16.6f}")
        txt.append(f"   Vacuum Volume: {self.vacuum_volume:16.6f}")
        txt.append(f"   No. of Electrons: {self.nelect:13.6f}")
        return '\n'.join(txt)

    @property
    def BCF(self):
        """
        Display the equivelent of the BCF.dat file.
        """
        idx = pd.Series(range(1,len(self.bader_df.index)+1))
        _df = self.bader_df.rename_axis(columns='#')
        _df = _df.set_index(idx)
        txt = _df.to_string(justify='center').split('\n')
        for i,t in enumerate(txt):
            txt[i] = ' ' + t
        txt.insert(1,'-'*(len(txt[1])+1))
        txt.append('-'*len(txt[1]))
        return '\n'.join(txt)

    def write(self,prefix=''):
        """
        Write ACF.dat and BCF.dat files.

        args:
            prefix (str): Prefix to ACF and BCF files.
        """
        with open(prefix+'ACF.dat','w') as f: f.write(self.ACF)
        with open(prefix+'BCF.dat','w') as f: f.write(self.BCF)
