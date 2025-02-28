from . import aux_functions as aux
import scipy as sc
import numpy as np
import random


def initialize_geometric_CPC(nx, ny, CPC_width=20, cohesin_width=4):
    phi = aux.dmatrix(nx, ny)
    for i in range(nx):
        for j in range(ny):
            if i > round(nx / 2) - CPC_width and i < round(nx / 2) + CPC_width:
                if j > round(ny / 2) - CPC_width and j < round(ny / 2) + CPC_width:
                    phi[i, j] = 1
                elif (
                    i > round(nx / 2) - cohesin_width
                    and i < round(nx / 2) + cohesin_width
                ):
                    phi[i, j] = 1
                else:
                    phi[i, j] = -1
            else:
                phi[i, j] = -1
    return phi


def initialize_round_CPC(nx, ny, CPC_width=10, cohesin_width=4):
    # Create an empty matrix filled with -1
    phi = np.ones(nx, ny)
    phi = -1 * phi

    # Define the center of the matrix
    center = (nx) / 2

    # Loop through each element of the matrix
    for i in range(nx):
        for j in range(ny):
            # Calculate the distance from the center
            distance = np.linalg.norm([i - center, j - center])

            # Check if the distance is less than or equal to CPC_width
            if distance <= CPC_width:
                phi[i, j] = 1.0
            elif (
                i > round((nx) / 2) - cohesin_width
                and i < round((nx) / 2) + cohesin_width
            ):
                phi[i, j] = 1.0
    return phi


def initialization_from_def(nx, ny, h, R0=0.1, epsilon=0.01):
    # Assuming aux.dmatrix(nx, ny) creates a zero matrix
    phi = np.zeros((ny, nx))

    x = np.arange(nx) * h
    y = np.arange(ny) * h

    # Manually creating the meshgrid effect
    R = np.sqrt((x[None, :] - 0.5) ** 2 + (y[:, None] - 0.5) ** 2)

    phi = np.tanh((R0 - R) / (np.sqrt(2) * epsilon))

    return phi


def initialization_random(nx, ny):
    return 2 * np.random.rand(nx, ny) - 1


def initialization_spinodal(nx, ny):
    return np.random.choice([-1, 1], size=(nx, ny))


def initialization_from_file(file, nx, ny, delim=",", transpose_matrix=False):
    phi = np.loadtxt(file, delim)
    if phi.shape != (nx, ny):
        print(
            "Warning: phi from file is wrong size: $(size(phi)) Expected: $(nx), $(ny)"
        )
    if transpose_matrix:
        phi = phi.transpose()

    return phi


def ch_initialization(
    nx,
    ny,
    method="spinodal",
    initial_file="",
    delim=",",
    h=1 / 128,
    R0=0.1,
    epsilon=0.01,
    cohesin_width=4,
    CPC_width=20,
):
    if method == "random":
        phi0 = initialization_random(nx, ny)
    elif method == "droplet":
        phi0 = initialization_from_def(nx, ny, h, R0=R0, epsilon=epsilon)
    elif method == "geometric":
        phi0 = initialize_geometric_CPC(
            nx, ny, CPC_width=CPC_width, cohesin_width=cohesin_width
        )
    elif method == "file":
        phi0 = initialization_from_file(initial_file, nx, ny, delim=delim)
    elif method == "spinodal":
        phi0 = initialization_spinodal(nx, ny)
    else:
        print(
            "Warning: initialize must be one of [random, droplet, geometric, file, spinodal]."
        )

    return phi0
