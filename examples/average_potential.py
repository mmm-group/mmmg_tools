import matplotlib.pyplot as plt
import numpy as np

from mmmg_tools.objects import Potential

pot = Potential.from_file("HARTREE.cube", ftype="cube")
z_av = pot.get_average_along_axis(2)
z_dist = np.linspace(0, pot.structure.lattice.c, pot.dim[2], False)

plt.plot(z_dist, z_av)
plt.xlabel(r"Distance ($\mathrm{\AA}$)")
plt.ylabel("Potential (eV)")
plt.show()
