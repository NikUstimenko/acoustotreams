import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0s = 2 * np.pi * np.linspace(1000, 300000, 100) / 343
materials = [(798, 1050), (1050, 2350), (998, 1497)]
thickness = 0.005
tr = np.zeros((2, len(k0s)))

def compute_coeffs(k0):
    pwb = acoustotreams.ScalarPlaneWaveBasisByComp.default([0, 0.1 * k0])
    slab = acoustotreams.AcousticSMatrices.slab(thickness, pwb, k0, materials)
    return slab.tr([1])
res = Parallel(n_jobs=-1)(delayed(compute_coeffs)(k0s[i]) for i in range(len(k0s)))
for i in range(len(k0s)):
    tr[0, i] = res[i][0][0]
    tr[1, i] = res[i][1][0]

fig, ax = plt.subplots()
ax.set_xlabel("Frequency (kHz)")
ax.plot(343 * k0s / (2 * np.pi) / 1000, tr[0, :])
ax.plot(343 * k0s / (2 * np.pi) / 1000, tr[1, :])
ax.legend(["$T$", "$R$"])
plt.show()