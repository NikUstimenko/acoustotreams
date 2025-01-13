import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0s =  2 * np.pi * np.linspace(10000, 50000, 200) / 343
material_slab = acoustotreams.AcousticMaterial(698, 950)
thickness = 0.0025
materials = [acoustotreams.AcousticMaterial(1050, 2350), 
            acoustotreams.AcousticMaterial(998, 1497)]
lmax = 4
radius = 0.0075
period = 0.025
lattice = acoustotreams.Lattice.square(period)


tr = np.zeros((len(k0s), 2))
def compute_coeffs(k0):
    kpar = [0, 0.1 * k0]
    spheres = acoustotreams.AcousticTMatrix.sphere(
        lmax, k0, radius, materials
        ).latticeinteraction.solve(lattice, kpar)
    
    pwb = acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders(kpar, lattice, 0.02)
    plw = acoustotreams.plane_wave_scalar(kpar, k0=k0, basis=pwb, material=materials[1])
    slab = acoustotreams.AcousticSMatrices.slab(thickness, pwb, k0, [materials[1], material_slab, materials[1]])
    dist = acoustotreams.AcousticSMatrices.propagation([0, 0, radius], pwb, k0, materials[1])
    array = acoustotreams.AcousticSMatrices.from_array(spheres, pwb)
    total = acoustotreams.AcousticSMatrices.stack([slab, dist, array])
    return total.tr(plw)

res = Parallel(n_jobs=-1)(delayed(compute_coeffs)(k0s[i]) for i in range(len(k0s)))
for i in range(len(k0s)):
    tr[i, 0] = res[i][0][0]
    tr[i, 1] = res[i][1][0]

fig, ax = plt.subplots()
ax.set_xlabel("Frequency (kHz)")
ax.plot(343 * k0s / (2 * np.pi) / 1000, tr[:, 0])
ax.plot(343 * k0s / (2 * np.pi) / 1000, tr[:, 1])
ax.legend(["$T$", "$R$"])
plt.show()