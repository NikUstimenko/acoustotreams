import acoutreams
import treams
import numpy as np

materials = [acoutreams.AcousticMaterial(1200, 1400), acoutreams.AcousticMaterial()]
inc = acoutreams.plane_wave_scalar([0, 0, 1], k0=1, material=acoutreams.AcousticMaterial())
tm = acoutreams.AcousticTMatrix.sphere(2, 1, 0.5, materials)
sca = tm.sca(inc)
#print(tm @ inc.expand(tm.basis, "regular"))
#print(tm.sca(inc))
#print(tm.sca(inc.expand(tm.basis, "regular")))

#basis = acoutreams.ScalarSphericalWaveBasis.default(2, positions=[0, 0, 1])
#inc = acoutreams.spherical_wave_scalar(1, 0, k0=1, basis=basis, material=acoutreams.AcousticMaterial(), modetype="regular")
#print(tm @ inc.expand(tm.basis, "regular"))
#print(tm.sca(inc))
#print(tm.sca(inc.expand(tm.basis, "regular")))

#grid = np.array([[0,0,1], [1, 0, 0]])
grid = np.array([0,0,1])

print(sca.pfield(grid))
puk

print(treams.lattice.lsumsw1d(1, 2, 4, 3, 0, 0))
acoutreams.ssw_Psi
acoutreams.Expand()

materials = [(1000, 1200), (1200, 1400)]
print()
print(acoutreams.lattice.lsumsw1d(1, 2, 4, 3, 0, 0))