import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

import acoustotreams

k0 =  2 * np.pi * 50000 / 343
materials = [acoustotreams.AcousticMaterial(1050 + 100j, 2350 - 300j), 
            acoustotreams.AcousticMaterial(998, 1497)]
lmax = 3
radii = [0.0075, 0.0065]
positions = [[-0.004, 0, -0.0075], [0.004, 0, 0.0075]]
period = 0.035
lattice = acoustotreams.Lattice(period)
kz = 0

spheres = [acoustotreams.AcousticTMatrix.sphere(lmax, k0, r, materials) for r in radii]
chain = acoustotreams.AcousticTMatrix.cluster(spheres, positions).latticeinteraction.solve(lattice, kz)

inc = acoustotreams.plane_wave_scalar([chain.k0, 0, 0], k0=chain.k0, material=chain.material)
sca = chain.sca(inc)

x = np.linspace(-0.75*period, 0.75*period, 101)
z = np.linspace(-0.5*period, 0.5*period, 101)
def compute_pressure(i, j):
    r = [x[j], 0, z[i]]  
    if chain.valid_points(r, radii):
        swb = acoustotreams.ScalarSphericalWaveBasis.default(0, positions=[r])
        result = sca.expandlattice(basis=swb).pfield(r) 
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

def compute_velocity(i, j):
    r = [x[j], 0, z[i]]  
    if chain.valid_points(r, radii):
        swb = acoustotreams.ScalarSphericalWaveBasis.default(1, positions=[r])
        result = sca.expandlattice(basis=swb).vfield(r)[0]  
    else:
        result = np.nan
    return i, j, result  
results = Parallel(n_jobs=-1)(
    delayed(compute_velocity)(i, j) 
    for i in range(len(z)) 
    for j in range(len(x))
)
vx = np.zeros((len(z), len(x)), complex)
for i, j, result in results:
    vx[i, j] = result

fig, ax = plt.subplots()
cax = ax.imshow(
    vx.real,
    extent = [x.min() * 100, x.max() * 100, z.min() * 100, z.max() * 100],
    aspect='equal',
    origin='lower',
)
cb = plt.colorbar(cax)
cb.set_label("Velocity $v_x$")
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