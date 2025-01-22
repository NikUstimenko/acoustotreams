import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0 =  2 * np.pi * 50000 / 343
materials = [acoustotreams.AcousticMaterial(1050 + 100j, 2350 - 300j), 
            acoustotreams.AcousticMaterial(998, 1497)]
lmax = mmax = 3
radii = [0.0075, 0.0065]
positions = [[-0.004, 0, -0.0075], [0.004, 0, 0.0075]]
period = 0.035
lattice = acoustotreams.Lattice(period)
kz = 0.1 * k0
n = 1497 / 343

spheres = [acoustotreams.AcousticTMatrix.sphere(lmax, k0, r, materials) for r in radii]
chain = acoustotreams.AcousticTMatrix.cluster(spheres, positions).latticeinteraction.solve(lattice, kz / n)

bmax = 3.1 * lattice.reciprocal
cwb = acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(kz / n, mmax, lattice, bmax, 2, positions)
chain_tmc = acoustotreams.AcousticTMatrixC.from_array(chain, cwb)

inc = acoustotreams.plane_wave_scalar(
    [np.sqrt(chain.k0 * chain.k0 - kz * kz), 0, kz], 
    k0=chain.k0, 
    material=chain.material
)
sca = chain.sca(inc)
sca_tmc = chain_tmc.sca(inc)

x = np.linspace(-0.75*period, 0.75*period, 101)
z = np.linspace(-0.5*period, 0.5*period, 101)
def compute_pressure(i, j):
    r = [x[j], 0, z[i]]  
    if chain_tmc.valid_points(r, radii):
        result = sca_tmc.pfield(r) 
    else:
        result = np.nan
        if chain.valid_points(r, radii):
            swb = acoustotreams.ScalarSphericalWaveBasis.default(0, positions=[r])
            result = sca.expandlattice(basis=swb).pfield(r) 
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
        np.real(p * np.exp(-1j * kz / n * period)),
        p.real,
        np.real(p * np.exp(1j * kz / n * period)),
    )
)
z = np.concatenate((z - period, z, z + period,))

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