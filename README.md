# Code for solving the Cahn Hilliard equation

Three options to solve the Cahn Hilliard equation using numerical solvers
1. Finite difference solver, adapted from https://github.com/nsbalbi/Spinodal-Decomposition. We found some problems in the original code, particularly in the laplacian function (which used a 9-point stencil and was rescaled incorrectly).
2. SAV solver
3. Multigrid (MG) solver


Either SAV or MG should be used for solving the CH equation, depending on context. See `run_spinodal_decomp.m` for an example of running all three methods.


To run this github repository:
1. Clone the repo to your local computer. (`git clone remote_address`)
2. Before making changes, always `git pull origin main` to get the latest version.
3. Make changes to the code.
4. `git add .` to add stage all changes for commit.
5. `git commit -m "commit message"` To commit your changes locally.
6. `git push origin main` to push your local changes to the remote repository. You will need to ensure there are no conflicts. This will open a pull request. If you are sure there are no conflicts, you can merge your new changes with the current `main` branch. 
