import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0s =  2 * np.pi * np.linspace(1000, 60000, 1000) / 343
material_slab = acoustotreams.AcousticMaterial(698, 950)
thickness = 0.0025
period = 0.025
lattice = acoustotreams.Lattice.square(period)
materials = [acoustotreams.AcousticMaterial(1050, 2350), 
            acoustotreams.AcousticMaterial(998, 1497)]
n = 1497 / 343
radius = 0.0075
lmax = mmax = 3

tr = np.zeros((len(k0s), 2))
def compute_coeffs(k0):
    kpar = np.array([0, 0.1 * k0])

    spheres = acoustotreams.AcousticTMatrix.sphere(
        lmax, k0, radius, materials
        ).latticeinteraction.solve(period, kpar[0] / n)
    
    bmax = 3.1 * 2 * np.pi / period
    cwb = acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(kpar[0] / n, mmax, period, bmax)
    spheres_tmc = acoustotreams.AcousticTMatrixC.from_array(
        spheres, cwb
    ).latticeinteraction.solve(period, kpar[1] / n)

    spwb = acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders(kpar / n, lattice, bmax)
    splw = acoustotreams.plane_wave_scalar(kpar / n, k0=k0, basis=spwb, material=materials[1])
    slab = acoustotreams.AcousticSMatrices.slab(thickness, spwb, k0, [materials[1], material_slab, materials[1]])
    dist = acoustotreams.AcousticSMatrices.propagation([0, 0, radius], spwb, k0, materials[1])
    array = acoustotreams.AcousticSMatrices.from_array(spheres_tmc, spwb)
    total = acoustotreams.AcousticSMatrices.stack([slab, dist, array])
    return total.tr(splw)

res = Parallel(n_jobs=-1)(delayed(compute_coeffs)(k0s[i]) for i in range(len(k0s)))
for i in range(len(k0s)):
    tr[i, 0] = res[i][0]
    tr[i, 1] = res[i][1]

fig, ax = plt.subplots()
ax.set_xlabel("Frequency (kHz)")
ax.plot(343 * k0s / (2 * np.pi) / 1000, tr[:, 0])
ax.plot(343 * k0s / (2 * np.pi) / 1000, tr[:, 1])
ax.legend(["$T$", "$R$"])
plt.show()