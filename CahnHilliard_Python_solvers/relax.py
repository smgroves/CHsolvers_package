import numpy as np
import CHsolvers

i = np.identity(3)  # identity matrix
mu = np.zeros((3, 3))
su = np.array([[118, -18, 0], [-18, 136, -18], [0, -18, 118]])
sw = np.zeros((3, 3))
c_relax = 2
xright = 0
xleft = 1
dt = 0.01
Cahn = 0.06**2
c_new_py, mu_new_py = CHsolvers.CahnHilliard_NMG(
    i, mu, su, sw, 3, 3, c_relax=c_relax, xright=xright, xleft=xleft, dt=dt, Cahn=Cahn
)
c_new_c = np.array(
    [
        [0.985498, 0.020633, -0.010602],
        [0.020633, 0.967676, 0.020633],
        [-0.010602, 0.020633, 0.985498],
    ]
)
mu_new_c = np.array(
    [
        [1.019644, -0.060925, -0.00202],
        [-0.060925, 1.028865, -0.060925],
        [-0.002025, -0.060925, 1.019644],
    ]
)

print("c_new_py = ", c_new_py)
print("c_new_c = ", c_new_c)
print("mu_new_py = ", mu_new_py)
print("mu_new_c = ", mu_new_c)
