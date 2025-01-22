import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0 =  2 * np.pi * 50000 / 343
materials = [acoustotreams.AcousticMaterial(1050 + 100j, 2350 - 300j), 
            acoustotreams.AcousticMaterial(998, 1497)]
lmax = mmax = 3
radii = [0.0055, 0.0045]
positions = [[-0.004, -0.005, 0], [0.004, 0.005, 0]]
lattice = 0.02
kz = 0.1 * k0
n = 1497 / 343

sphere = acoustotreams.AcousticTMatrix.sphere(lmax, k0, radii[0], materials)
chain = sphere.latticeinteraction.solve(lattice, kz / n)
bmax = 3.1 * 2 * np.pi / lattice
cwb = acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(kz / n, mmax, lattice, bmax)
chain_tmc = acoustotreams.AcousticTMatrixC.from_array(chain, cwb)

cylinder = acoustotreams.AcousticTMatrixC.cylinder(np.unique(cwb.kz), mmax, k0, radii[1], materials)

cluster = acoustotreams.AcousticTMatrixC.cluster([chain_tmc, cylinder], positions).interaction.solve()
inc = acoustotreams.plane_wave_scalar(
    [np.sqrt(cluster.k0 * cluster.k0 - kz * kz), 0, kz], 
    k0=cluster.k0, 
    material=cluster.material
)
sca = cluster.sca(inc)

x = np.linspace(-lattice, lattice, 101)
z = np.linspace(-0.5*lattice, 0.5*lattice, 101)
def compute_pressure(i, j):
    r = [x[j], 0, z[i]]  
    if cluster.valid_points(r, radii):
        result = sca.pfield(r)
        print(result) 
    else:
        result = np.nan
    return i, j, result  
results = Parallel(n_jobs=-1)(
    delayed(compute_pressure)(i, j) 
    for i in range(len(z)) 
    for j in range(len(x))
)
p = np.zeros((len(z), len(x)), complex)
for i, j, result in results:
    p[i, j] = result

p = np.concatenate(
    (
        np.real(p * np.exp(-1j * kz / n * lattice)),
        p.real,
        np.real(p * np.exp(1j * kz / n * lattice)),
    )
)
z = np.concatenate((z - lattice, z, z + lattice,))

fig, ax = plt.subplots()
cax = ax.imshow(
    p,
    extent = [x.min() * 100, x.max() * 100, z.min() * 100, z.max() * 100],
    aspect='equal',
    origin='lower',
)
cb = plt.colorbar(cax)
cb.set_label("Pressure")
ax.set_xlabel("x (cm)")
ax.set_ylabel("z (cm)")
ax.annotate(
    "", 
    xy=(-1.0, 0),
    xytext=(-1.75, 0), 
    arrowprops=dict(
        arrowstyle="->",
        lw=3, 
        color="red" 
    )
)
plt.show()