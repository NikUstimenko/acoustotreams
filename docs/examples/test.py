import numpy as np
import acoustotreams
import treams

print(acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, 0], acoustotreams.Lattice.square(2 * np.pi), 1))