import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0 =  2 * np.pi * 50000 / 343
materials = [acoustotreams.AcousticMaterial(1050 + 100j, 2350 - 300j), 
            acoustotreams.AcousticMaterial(998, 1497)]
lmax = mmax = 3
radii = [0.0065, 0.0055]
positions = [[-0.004, -0.005, 0], [0.004, 0.005, 0]]
period = 0.035
lattice_z = acoustotreams.Lattice(period)
lattice_x = acoustotreams.Lattice(period, "x")
kz = 0.1 * k0
n = 1497 / 343

sphere = acoustotreams.AcousticTMatrix.sphere(lmax, k0, radii[0], materials)
chain = sphere.latticeinteraction.solve(lattice_z, kz / n)
bmax = 3.1 * lattice_z.reciprocal
cwb = acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(kz / n, mmax, lattice_z, bmax)
kzs = np.unique(cwb.kz)
chain_tmc = acoustotreams.AcousticTMatrixC.from_array(chain, cwb)
cylinder = acoustotreams.AcousticTMatrixC.cylinder(kzs, mmax, k0, radii[1], materials)

cluster = acoustotreams.AcousticTMatrixC.cluster(
    [chain_tmc, cylinder], positions
).latticeinteraction.solve(lattice_x, 0)
inc = acoustotreams.plane_wave_scalar(
    [0, np.sqrt(cluster.k0 * cluster.k0 - kz * kz), kz], 
    k0=cluster.k0, 
    material=cluster.material
)
sca = cluster.sca(inc)

x = np.linspace(-0.5*period, 0.5*period, 101)
z = np.linspace(-0.5*period, 0.5*period, 101)
def compute_pressure(i, j):
    r = [x[j], 0, z[i]]  
    if cluster.valid_points(r, radii):
        cwb = acoustotreams.ScalarCylindricalWaveBasis.default(kzs, 0, positions=[r])
        result = sca.expandlattice(basis=cwb).pfield(r)
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

fig, ax = plt.subplots()
cax = ax.imshow(
    p.real,
    extent = [x.min() * 100, x.max() * 100, z.min() * 100, z.max() * 100],
    aspect='equal',
    origin='lower',
)
cb = plt.colorbar(cax)
cb.set_label("Pressure")
ax.set_xlabel("x (cm)")
ax.set_ylabel("z (cm)")
plt.show()