import numpy as np
from . import aux_functions_SAV as aux   # Import auxiliary functions

# This function uses the sav method to solve the
# Cahn-Hilliard equation for the next time step.

# INPUTS
# phi_old  = Prior chemical state.
# r_old    = Prior SAV parameter.
# hx, hy   = Grid spacings in x and y directions.
# k2, k4   = Laplacian and biharmonic operators.
# dt       = Time step.
# epsilon2 = Square of the interface width parameter.
# boundary = 'periodic' or 'neumann' (not explicitly used in this snippet).
# C0       = Regularization parameter.
# Beta     = Relaxation parameter.
# gamma0   = Stabilization parameter.

# OUTPUT
# phi_new = Next chemical state.
# r_new   = Next sav state.


def sav_solver(phi_old, r_old, hx, hy, k2, k4, dt, epsilon2, boundary, C0, Beta, gamma0):

    phi0 = phi_old
    r0 = r_old

    phi0_df = aux.df(phi0, gamma0)  # df at phi0
    Lap_dfphi0 = aux.Lap_SAV(phi0_df, k2)  # Lap of df(phi0)
    phi_bar = aux.A_inv_CN(phi0 + dt/2 * Lap_dfphi0,
                           dt, k2, k4, gamma0, epsilon2)

    # Step 1
    b = aux.b_fun(phi_bar, hx, hy, C0, gamma0)
    g = aux.g_fun_CN(phi0, r0, b, dt, hx, hy, epsilon2, gamma0, Beta, C0, k2)

    AiLb = aux.A_inv_CN(aux.Lap_SAV(b, k2), dt, k2, k4, gamma0, epsilon2)
    Aig = aux.A_inv_CN(g, dt, k2, k4, gamma0, epsilon2)

    gamma = -np.real(np.fft.fft2(b*AiLb))
    gamma = gamma[0, 0]*hx*hy

    # Step 2
    bphi = np.real(np.fft.fft2(b*Aig))
    bphi = bphi[0, 0]*hx*hy/(1+dt/4*gamma)

    # Step 3
    phi_new = dt/4*bphi*AiLb + Aig
    r_new = aux.r_fun(phi_new, phi_old, r0, b, hx, hy, C0, Beta, dt)
    return phi_new, r_new
