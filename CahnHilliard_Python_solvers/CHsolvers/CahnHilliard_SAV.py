import numpy as np
import os
from .aux_functions import *
from .aux_functions_SAV import *
from . import SAV_solver as solver


@time_and_mem
def CahnHilliard_SAV(phi0, t_iter=1e3, dt=2.5e-5, dt_out=1, m=4, epsilon2=np.nan, boundary="periodic", domain=[1, 0, 1, 0], printres=False, printphi=False, pathname="cd", C0=0, Beta=0, gamma0=0):
    """
    This function uses the nonlinear multigrid method to solve the Cahn-Hilliard equation for a specified number of time steps of size dt.

    REQUIRED INPUT
        phi0 = Initial field of chemical states in the domain, created by ch_initialization.

    OPTIONAL INPUT
        t_iter = Number of time steps simulated. Default = 1e3.
        dt = Time step. Default = 2.5e-5 characteristic times.
        solver_iter = Number of solver iterations per time step. Default = 1e4.
        tol = Solver tolerance per time step. Default = 1e-5.
        dt_out = Spacing of time steps output to phi_t as a multidimensional
           array (if less than 1e9 elements) or printed file (if greater than
           1e9 elements). Default = 1, to output every time step.
        m = Number of mesh points over which the interface exists. Default = 4.
        epsilon2 = Squared interfacial transition distance; if specified,
           m will be overwritten. Default = nan (do not overwrite m).
        boundary = Boundary conditions for the simulation:
           'periodic' (default) - flux on one domain border equals negative flux on the opposite border.
           'neumann' - zero flux on the domain borders.
        c_relax = Number of smoothing relaxations done at the start and
           of each V-cycle. Default = 2.
        domain = Vector of rightmost and leftmost grid points in x and y.
           Format: [xright xleft yright yleft]. Default = [1 0 1 0].
        printres = Logical to print solver residuals to a file. Default = False.
        printphi = Logical to print phi to a file regardless of whether or
           not it can be saved as a multidimensional array. Default = False.
         pathname = Name of the path to which phi is printed. Default = 'cd'. May include a prefix for the filename.

    OUTPUT
        t_out = Time corresponding to the dt time step outputs.
        phi_t = Multidimensional array of phi over t_out.
        delta_mass_t = Vector of mass change over t_out.
        E_t = Vector of relative energy over t_out.
    """
    [nx, ny] = phi0.shape
    xright, xleft, yright, yleft = domain
    Lx = xright - xleft
    Ly = yright - yleft

    # % Decide on the solver's mesh spacing for NEUMANN vs PERIODIC
    # %  - For Neumann: we will mirror the domain, so pass 2*hx and 2*hy into sav_solver.
    # %  - For Periodic: keep as-is.
    if boundary == "neumann":
        Lx = 2 * Lx
        Ly = 2 * Ly
        nx = 2 * nx
        ny = 2 * ny
    hx = Lx / nx
    hy = Ly / ny
    h2 = hx * hy  # Define mesh size

    if np.isnan(epsilon2):
        # Define ϵ^2 if not prespecified
        epsilon2 = h2 * m ** 2 / (2 * np.sqrt(2) * np.arctanh(0.9)) ** 2
    else:
        # Else overwrite m
        m = np.sqrt((epsilon2 * (2 * np.sqrt(2) * np.arctanh(0.9)) ** 2) / h2)

    kx = compute_kx(nx, Lx)
    ky = compute_kx(ny, Ly)

    kxx = kx**2
    kyy = ky**2
    kxx_mat, kyy_mat = meshgrid(kxx, kyy)

    k2 = kxx_mat + kyy_mat
    k4 = k2**2

    if boundary == "neumann":
        phi_old = ext(phi0)
    elif boundary == "periodic":
        phi_old = phi0.copy()

    r_old = r0_fun(phi_old, hx, hy, C0)  # Initialize sav state

    # Logical index for the need to downsample
    downsampled = nx * ny * t_iter / dt_out > 1e9
    n_timesteps = np.floor(t_iter / dt_out)

    if printphi:
        mass_t = np.zeros(int(n_timesteps + 1))
        E_t = np.zeros(int(n_timesteps + 1))
        t_out = np.arange(0, (int(n_timesteps+1))*dt*dt_out, dt_out*dt)
        if pathname == "cd":
            pathname = os.pwd()

        if boundary == "neumann":
            phi_old_out = extback(phi_old)
        else:
            phi_old_out = phi_old

        with open(f"{pathname}phi.csv", "w") as f:
            for i in range(nx):
                for j in range(ny):
                    f.write(f"{phi_old_out[i][j]},")
                f.write("\n")  # end of each row

        phi_t = np.nan
    else:  # save to variables
        if downsampled:
            # we need to round up to ensure we have enough space
            new_dt_out = np.ceil(nx * ny * t_iter / 1e9)
            print(
                "Variable phi_t is too large with dt_out = %4.0f. Downsampling to every %4.0f time steps\n", dt_out, new_dt_out)
            dt_out = new_dt_out
            n_timesteps = np.floor(t_iter / dt_out)

        mass_t = np.zeros(int(n_timesteps + 1))
        E_t = np.zeros(int(n_timesteps + 1))
        t_out = np.arange(0, (int(n_timesteps+1))*dt*dt_out, dt_out*dt)

        if boundary == "neumann":
            phi_t = np.zeros((int(nx / 2), int(ny / 2), int(n_timesteps + 1)))
            phi_old_out = extback(phi_old)
            phi_t[:, :, 1] = phi_old_out
        else:
            phi_t = np.zeros((nx, ny, int(n_timesteps + 1)))
            phi_old_out = phi_old
            phi_t[:, :, 1] = phi_old_out

    mass_t[1] = calculate_mass(phi0, h2, nx, ny)
    E_t[1] = calculate_discrete_energy(phi0, h2, nx, ny, epsilon2)
    if printres:
        print(
            "Saving squared residuals per iteration to file in the output directory\n")

    for it in range(t_iter):
        phi_new, r_new = solver.sav_solver(
            phi_old, r_old, hx, hy, k2, k4, dt, epsilon2, boundary, C0, Beta, gamma0)
        if boundary == "neumann":
            phi_new_out = extback(phi_new)
        else:
            phi_new_out = phi_new

        if it / t_iter * 100 % 5 == 0:
            print((f"{it / t_iter * 100:3.0f} percent complete"))

        if it % dt_out == 0:
            t_index = int(np.floor(it / dt_out) + 1)
            if printphi:
                with open(f"{pathname}phi.csv", "a") as f:
                    for i in range(nx):
                        for j in range(ny):
                            f.write(f"{phi_new_out[i][j]},")
                        f.write("\n")  # end of each row

            else:
                phi_t[:, :, t_index] = phi_new_out

            mass_t[t_index] = calculate_mass(phi_new_out, h2, nx, ny)
            E_t[t_index] = calculate_discrete_energy(
                phi_new_out, h2, nx, ny, epsilon2)

        phi_old = (phi_new).copy()
        r_old = (r_new).copy()

    delta_mass_t = mass_t - mass_t[1]
    E_t = E_t / E_t[1]
    return t_out, phi_t, delta_mass_t, E_t
