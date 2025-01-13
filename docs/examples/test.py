import numpy as np
import acoustotreams
import treams

plw = acoustotreams.plane_wave_scalar([0, 3, 4], k0=5, material=[(1000, 2000, 3000)])
print(plw.material)