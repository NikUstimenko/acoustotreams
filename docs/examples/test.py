import numpy as np
import acoutreams
import treams

sscw = acoutreams.cylindrical_wave_scalar(0, 1, k0=1, material=1.3)
ex = acoutreams.ExpandLattice(basis=acoutreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, .1], acoutreams.Lattice([7, 7], "zx"), 1))
print((ex @ sscw).basis)