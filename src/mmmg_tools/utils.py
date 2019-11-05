import numpy as np

class sparray():
    """
    Sparse array class for nd-arrays.
    """
    def __init__(self,shape,data=None):
        """
        Initialise the object. 
        
        Data passed must be 2D array of shape (2,N) where (0,:) is an 
        array of int corresponding to the index of the flattened array 
        and (1,:) is the value at each index.

        args:
            shape (3x1 array): Dimensions of ndarray.
            data (2xN array): The data of the sparse array.
        """
        self.shape = shape
        self._data = {}
        if data is not None:
            self._data = {int(k): v for k, v in data.T}

    def __repr__(self):
        return f'{self.__class__}({self.__dict__})'

    def __setitem__(self,idx,v):
        if idx is tuple:
            idx = np.ravel_multi_index(idx,self.shape)
        self._data[idx] = v

    def __getitem__(self,idx):
        if idx is tuple:
            idx = np.ravel_multi_index(idx,self.shape)
        return self._data[idx]

    def __delitem__(self,idx):
        if idx is tuple:
            idx = np.ravel_multi_index(idx,self.shape)
        del(self._data[idx])

    @property
    def sum(self):
        """
        Returns the sum of the sparse array.
        """
        c = 0
        for v in self._data.values():
            c += v
        return c

    @property 
    def pop_percent(self):
        """
        Returns the population percentage of the sparse array.
        """
        return f'{(100 * self.sum) / np.prod(self.shape):4.4} %'

    def shape_match(self, other):
        """
        Check shape of arrays.
        """
        err = f'Array shapes do not match: {self.shape} {other.shape}'
        if self.shape != other.shape:
            raise ValueError(err)

    def copy(self):
        """
        Copy sparse array.
        """
        out = self.__class__(self.shape)
        out._data = self._data.copy()
        return out

    def __add__(self, other):
        """
        Adds an array to a sparse array.
        """
        self.shape_match(other)
        out = other.copy()
        if other.__class__ == self.__class__:
            for k,v in self._data.items():
                out[k] = other._data.get(k,0) + v
        else:
            for k,v in self._data.items():
                idx = np.unravel_index(k,self.shape)
                out[idx] += v
        return out

    def __iadd__(self, other):
        """
        Adds an array to a sparse array.
        """
        self.shape_match(other)
        out = other.copy()
        if other.__class__ == self.__class__:
            for k,v in self._data.items():
                out[k] = other._data.get(k,0) + v
        else:
            for k,v in self._data.items():
                idx = np.unravel_index(k,self.shape)
                out[idx] += v
        return out

    def __sub__(self, other):
        """
        Subtracts an array from a sparse array.
        """
        self.shape_match(other)
        out = other.copy()
        if other.__class__ == self.__class__:
            for k,v in other._data.items():
                out[k] = -v
        else:
            out *= -1
        return self.__add__(out)

    def __mul__(self, other):
        """
        Multiplies a sparse array with an array. Always returns a 
        sparse array.
        """
        self.shape_match(other)
        out = self.copy()
        if other.__class__ == self.__class__:
            for k,v in self._data.items():
                if other._data.has_key(k):
                    out[k] = other[k] * v
                else:
                    del(out[k])
        else:
            for k,v in self._data.items():
                idx = np.unravel_index(k,self.shape)
                if other[idx] != 0:
                    out[k] = v * other[idx]
        return out

    def __div__(self, other):
        """
        Divides a sparse array by an array. Always returns a sparse
        array.
        """
        self.shape_match(other)
        out = self.__class__(self.shape)
        if self.__class__ == other.__class__:
            for k,v in other._data.items():
                    out[k] = 1/v
        else:
            for k,_ in self._data.items():
                idx = np.unravel_index(k,self.shape)
                if other[idx] != 0:
                    out[k] = 1/other[idx]
        return self.__mul__(out)

    def fill(self):
        """
        Create dense numpy array from sparse array.
        """
        out = np.zeros(self.shape)
        return self.__add__(out)

    def apply(self, obj):
        """
        Apply the mask to a Charge or Potential object returning new object.

        args:
            obj (object): Charge or Potential object. (Doesn't accept pymatgen objects)
        """
        data = {k: self.__mul__(v).fill() for k, v in obj.data.items()}
        return obj.__class__(obj.to_structure(), data)

def rotate_to_vector(vector:list):
    """
    Calculate 3D rotation matrix for aligning z-vector and a vector.

    args:
        vector: 3D vector != [0, 0, 1]

    returns:
        Q: rotation matrix
    """
    r = np.array(vector) / np.linalg.norm(vector)
    theta, b = np.arccos(r[2]), np.cross([0,0,1], r)
    b = b / np.linalg.norm(b)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    q = [c,*np.multiply(s, b)]
    Q = np.array([
        [
            q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 
            2 * (q[1] * q[2] - q[0] * q[3]), 
            2 * (q[1] * q[3] + q[0] * q[2]),
        ],
        [
            2 * (q[1] * q[2] + q[0] * q[3]), 
            q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2, 
            2 * (q[2] * q[3] - q[0] * q[1]),
        ],
        [
            2 * (q[1] * q[3] - q[0] * q[2]),
            2 * (q[2] * q[3] + q[0] * q[1]),
            q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2,
        ],
    ])
    return Q
