import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0s = 2 * np.pi * np.linspace(1000, 115000, 400) / 343
materials = [acoustotreams.AcousticMaterial(1050, 2350), acoustotreams.AcousticMaterial(998, 1497)]
mmax = 6
radius = 0.005
lattice = acoustotreams.Lattice.cubic(0.015)
kpar = [0, 0, 0]

def compute_svd(k0):
    cyl = acoustotreams.AcousticTMatrixC.cylinder(kpar[2], mmax, k0, radius, materials)
    svd = np.linalg.svd(cyl.latticeinteraction(lattice, kpar[:2]), compute_uv=False)
    return svd[-1]
res = Parallel(n_jobs=-1)(delayed(compute_svd)(k0s[i]) for i in range(len(k0s)))

fig, ax = plt.subplots()
ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("Smallest singular value")
ax.semilogy(343 * k0s / (2 * np.pi) / 1000, res)
fig.show()