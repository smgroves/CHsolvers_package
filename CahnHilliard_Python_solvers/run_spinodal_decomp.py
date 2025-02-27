import CHsolvers as ch
import pandas as pd
import os
import time
ch.init.initialization(100, 100)
start = time.time()
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
