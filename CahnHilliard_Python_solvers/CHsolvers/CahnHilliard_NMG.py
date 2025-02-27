import numpy as np
import random
import os
import time
import pandas as pd
from . import NMG_solver as solver
from . import initialization as init

######################
# Global variables
######################
nx = 64  # max number of grid points in x-direction defined as a global variable
ny = 64  # max number of grid points in y-direction defined as a global variable

n_level = int(np.log(nx) / np.log(2.0) + 0.1)  # original c code uses natural log too
c_relax = int(2)  # number of SMOOTH relaxation operations defined as a global variable
xleft = 0.0  # left x-coordinate defined as a global variable
xright = 1.0  # right x-coordinate defined as a global variable
yleft = 0.0  # left y-coordinate defined as a global variable
yright = 1.0  # right y-coordinate defined as a global variable

# todo: check if these are needed; it appears that only ht2 (temp h^2) is used in the code
h = xright / float(nx)  # space step size defined as a global variable
h2 = h**2  # space step size squared defined as a global variable
dt = 0.1 * h2  # ∆t defined as a global variable
gam = 4 * h / (2 * np.sqrt(2) * np.arctanh(0.9))
Cahn = gam**2  # ϵ^2 defined as a global variable


if __name__ == "__main__":
    for max_it in [10000]:
        for solver_iter in [1000, 10000, 100000]:
            brcd = random.randint(0, 1000)
            start = time.time()
            # solver_iter = 1
            tol = 1e-6
            # max_it = 100 # number of time steps
            ns = 10  # frequency that time step is printed to file (1 = every time step, 10 = every 10 steps, etc)
            print(
                f"nx = {nx}, ny = {ny}, dt = {dt}, Cahn = {Cahn}, max_it = {max_it}, ns = {ns}, n_level = {n_level}"
            )
            mu = np.zeros((nx, ny))  # µ defined as a global variable
            oc = init.initialization(nx, ny)
            # print(oc)
            nc = oc.copy()

            for it in range(max_it):
                nc = solver.cahn(
                    oc, nc, mu, solver_iter=solver_iter, tol=tol
                )  # update phi
                oc = nc.copy()  # copy updated phi to old phi

            with open(f"../outputs/runtime_tests/phi_{brcd}.txt", "w") as f:
                # write initial state to file
                for i in range(nx):
                    for j in range(ny):
                        f.write(f"{oc[i][j]} ")
                    f.write("\n")
                for it in range(max_it):
                    nc = solver.cahn(
                        oc, nc, mu, solver_iter=solver_iter, tol=tol
                    )  # update phi
                    oc = nc.copy()  # copy updated phi to old phi
                    if it % ns == 0:
                        # write state to file every ns iterations
                        for i in range(nx):
                            for j in range(ny):
                                f.write(f"{oc[i][j]} ")
                            f.write("\n")
                        print(it)
            end = time.time()
            print(f"Time elapsed: {end - start} seconds")
            T = {}
            # T['brcd'] = brcd
            T["c_relax"] = c_relax
            T["Cahn"] = Cahn
            T["solver"] = "Python"
            T["dt"] = dt
            T["max_it"] = max_it
            T["solver_iter"] = solver_iter
            T["n_level"] = n_level
            T["ns"] = ns
            T["nx"] = nx
            T["ny"] = ny
            T["time"] = end - start
            T["tol"] = tol
            T = pd.DataFrame([T])
            file = "../Job_specs_all_py_c.csv"
            if max_it != 1:
                if not os.path.isfile(file):
                    T.to_csv(file)
                else:
                    with open(file, "a") as f:
                        T.to_csv(f, header=False, index=False)
