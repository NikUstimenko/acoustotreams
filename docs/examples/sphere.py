import matplotlib.pyplot as plt
import numpy as np

import acoustotreams

k0s = 2 * np.pi * np.linspace(50000, 500000, 200) / acoustotreams.AcousticMaterial().c
materials = [acoustotreams.AcousticMaterial(1050 + 100j, 2350 - 300j), 
             acoustotreams.AcousticMaterial(998, 1497)]
lmax = 10
radius = 0.005
spheres = [acoustotreams.AcousticTMatrix.sphere(lmax, k0, radius, materials) for k0 in k0s]

xs_sca = np.array([tm.xs_sca_avg for tm in spheres]) / (np.pi * radius ** 2)
xs_ext = np.array([tm.xs_ext_avg for tm in spheres]) / (np.pi * radius ** 2)

swb_lmax1 = acoustotreams.ScalarSphericalWaveBasis.default(1)
spheres_lmax1 = [tm[swb_lmax1] for tm in spheres]
xs_sca_lmax1 = np.array([tm.xs_sca_avg for tm in spheres_lmax1]) / (np.pi * radius ** 2)
xs_ext_lmax1 = np.array([tm.xs_ext_avg for tm in spheres_lmax1]) / (np.pi * radius ** 2)

tm = spheres[50]
inc = acoustotreams.plane_wave_scalar([0, 0, tm.k0], k0=tm.k0, material=tm.material)
sca = tm.sca(inc)
x = np.linspace(-0.0075, 0.0075, 101)
z = np.linspace(-0.0075, 0.0075, 101)
X, Z = np.meshgrid(x, z)
points = np.stack([X, np.zeros_like(X), Z], axis=-1)
valid = tm.valid_points(points, [radius])
intensity = np.full(X.shape, np.nan)
field = inc.pfield(points[valid]) + sca.pfield(points[valid])
intensity[valid] = np.abs(np.asarray(field)[:, 0]) ** 2

theta = np.linspace(0, 2 * np.pi, 301)
directions = np.stack([np.sin(theta), np.zeros_like(theta), np.cos(theta)], axis=-1)
radpattern = np.abs(np.asarray(sca.pamplitudeff(directions))[:, 0]) ** 2

import matplotlib
matplotlib.use("Agg")

fig, ax = plt.subplots()
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xs_ext)
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xs_sca)
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xs_ext_lmax1, color="C0", linestyle=":")
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xs_sca_lmax1, color="C1", linestyle=":")
ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("Efficiency")
ax.legend(["Extinction", "Scattering", "Extinction $l_{\\rm max}=1$", "Scattering $l_{\\rm max}=1$"])
fig.show()

fig, ax = plt.subplots()
cax = ax.imshow(
    intensity,
    extent = [x.min() * 100, x.max() * 100, z.min() * 100, z.max() * 100],
    aspect='auto',
    origin='lower',
)
cb = plt.colorbar(cax)
cb.set_label(r"$|p_{\rm tot}|^2/|p_0|^2$")
ax.set_xlabel("x (cm)")
ax.set_ylabel("z (cm)")
ax.annotate(
    "", 
    xy=(0, -0.5),
    xytext=(0, -0.75), 
    arrowprops=dict(
        arrowstyle="->",
        lw=3, 
        color="red" 
    )
)
fig.show()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, radpattern / np.max(radpattern))
ax.set_theta_offset(np.pi / 2)
ax.annotate(
    "", 
    xy=(np.pi, 0.2), 
    xytext=(np.pi, 0.6),
    arrowprops=dict(
        arrowstyle="->", 
        color="red", 
        lw=3
    )
)
fig.show()