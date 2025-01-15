import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0s =  2 * np.pi * np.linspace(100, 100000, 800) / 343
material_slab = acoustotreams.AcousticMaterial(698, 950)
thickness = 0.0025
period = 0.025
lattice = acoustotreams.Lattice.square(period)
materials = [acoustotreams.AcousticMaterial(1050, 2350), 
            acoustotreams.AcousticMaterial(998, 1497)]
n = materials[1].c / 343 
radius = 0.0075
lmax = 4
az = 0.02

def compute_ndiff(period, wavelength):
    if period < wavelength:
        return 0.
    if period >= wavelength and period < np.sqrt(2) * wavelength:
        return 1.
    if period >= np.sqrt(2) * wavelength and period < 2 * wavelength:
        return np.sqrt(2)

def compute_coeffs(k0):
    kpar = np.array([0, 0])

    spheres = acoustotreams.AcousticTMatrix.sphere(
        lmax, k0, radius, materials
        ).latticeinteraction.solve(lattice, kpar / n)
    
    bmax = 1.1 * 2 * np.pi / period * compute_ndiff(period, 2 * np.pi * n / k0)
    spwb = acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders(kpar / n, lattice, bmax)
    slab = acoustotreams.AcousticSMatrices.slab(thickness, spwb, k0, [materials[1], material_slab, materials[1]])
    dist = acoustotreams.AcousticSMatrices.propagation([0, 0, radius], spwb, k0, materials[1])
    array = acoustotreams.AcousticSMatrices.from_array(spheres, spwb)
    total = acoustotreams.AcousticSMatrices.stack([slab, dist, array])
    x, _ = total.bands_kz(az)
    x = x * az / np.pi
    return x[np.abs(np.imag(x)) < 0.1]

res = Parallel(n_jobs=-1)(delayed(compute_coeffs)(k0s[i]) for i in range(len(k0s)))

fig, ax = plt.subplots()
for k0, sel in zip(k0s, res):
    ax.scatter(sel.real, len(sel) * [343 * k0 / (2 * np.pi) / 1000], 0.2, c="C0")
    ax.scatter(sel.imag, len(sel) * [343 * k0 / (2 * np.pi) / 1000], 0.2, c="C1")
ax.set_xlabel("$k_z a_z / \\pi$")
ax.set_ylabel("Frequency (kHz)")
ax.set_xlim([-1, 1])
ax.set_ylim(ymin=0)
ax.legend(["$Real$", "$Imag$"])
plt.show()