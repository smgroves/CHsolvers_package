import profile
import CHsolvers as ch
import numpy as np
import cProfile
print(cProfile.__file__)
phi0 = np.random.uniform(-0.1, 0.1, (64, 64))
profile.run('ch.NMG.CahnHilliard_NMG(phi0, t_iter=5)', sort='tottime')
