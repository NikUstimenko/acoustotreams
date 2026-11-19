
import matplotlib.pyplot as plt
import numpy as np
 
import acoustotreams
 
k0s = 2 * np.pi * np.linspace(50000, 500000, 200) / acoustotreams.AcousticMaterial().c
materials = [acoustotreams.AcousticMaterial(1050 + 100j, 2350 - 300j), 
             acoustotreams.AcousticMaterial(998, 1497)]
kzs = [0]
mmax = 10
radius = 0.005
cylinders = [acoustotreams.AcousticTMatrixC.cylinder(kzs, mmax, k0, radius, materials) for k0 in k0s]
 
xw_sca = np.array([tm.xw_sca_avg for tm in cylinders]) / (2 * radius)
xw_ext = np.array([tm.xw_ext_avg for tm in cylinders]) / (2 * radius)
 
cwb_mmax1 = acoustotreams.ScalarCylindricalWaveBasis.default(kzs, 1)
cylinders_mmax0 = [tm[cwb_mmax1] for tm in cylinders]
xw_sca_mmax0 = np.array([tm.xw_sca_avg for tm in cylinders_mmax0]) / (2 * radius)
xw_ext_mmax0 = np.array([tm.xw_ext_avg for tm in cylinders_mmax0]) / (2 * radius)
 
tm = cylinders[50]
inc = acoustotreams.plane_wave_scalar([tm.k0, 0, 0], k0=tm.k0, material=tm.material)
sca = tm.sca(inc)
 
x = np.linspace(-0.0075, 0.0075, 101)
y = np.linspace(-0.0075, 0.0075, 101)
X, Y = np.meshgrid(x, y)
points = np.stack([X, Y, np.zeros_like(X)], axis=-1)
valid = tm.valid_points(points, [radius])
intensity = np.full(X.shape, np.nan)
field = inc.pfield(points[valid]) + sca.pfield(points[valid])
intensity[valid] = np.abs(np.asarray(field)[:, 0]) ** 2
 
phi = np.linspace(0, 2 * np.pi, 301)
directions = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
radpattern = np.abs(np.asarray(sca.pamplitudeff(directions))[:, 0]) ** 2
 
fig, ax = plt.subplots()
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xw_ext)
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xw_sca)
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xw_ext_mmax0, color="C0", linestyle=":")
ax.plot(k0s * acoustotreams.AcousticMaterial().c / (2 * np.pi) / 1000, xw_sca_mmax0, color="C1", linestyle=":")
ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("Efficiency")
ax.legend(["Extinction", "Scattering", "Extinction $m_{\\rm max}=1$", "Scattering $m_{\\rm max}=1$"])
fig.show()
 
fig, ax = plt.subplots()
cax = ax.imshow(
    intensity,
    extent = [x.min() * 100, x.max() * 100, y.min() * 100, y.max() * 100],
    aspect='auto',
    origin='lower',
)
cb = plt.colorbar(cax)
cb.set_label("Intensity")
ax.set_xlabel("x (cm)")
ax.set_ylabel("y (cm)")
ax.annotate(
    "", 
    xy=(-0.5, 0),
    xytext=(-0.75, 0), 
    arrowprops=dict(
        arrowstyle="->",
        lw=3, 
        color="red" 
    )
)
fig.show()
 
ax = plt.figure()
ax = plt.subplot(111, polar=True) 
ax.plot(phi, radpattern / np.max(radpattern) )
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