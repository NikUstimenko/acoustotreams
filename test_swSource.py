import numpy as np
import treams

k0 = 2 * np.pi / 1000
materials = [treams.Material(16 + 0.5j), treams.Material()]
lmax = 2
radii = [110, 90, 80, 75]
positions = (220 / np.sqrt(24)) * np.array(
    [
        [np.sqrt(8), 0, -1],
        [-np.sqrt(2), np.sqrt(6), -1],
        [-np.sqrt(2), -np.sqrt(6), -1],
        [0, 0, 3],
    ]
)

spheres = [treams.TMatrix.sphere(lmax, k0, r, materials) for r in radii]
tm = treams.TMatrix.cluster(spheres, positions).interaction.solve()
off_centered_swb = treams.SphericalWaveBasis.default(2, positions=[[0, 0, -100]])
inc = treams.spherical_wave(1, 0, 0, k0=tm.k0, basis=off_centered_swb, material=tm.material, modetype="singular")
print(tm.xs(inc.expand(tm.basis, "regular")))
