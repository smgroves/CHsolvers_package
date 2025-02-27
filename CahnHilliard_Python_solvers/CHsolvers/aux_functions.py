import numpy as np


def dmatrix(nrows, ncols):
    """
    Create an empty matrix of size nrows x ncols
    :param nrows: Number of rows
    :param ncols: Number of columns
    :return: Matrix of size nrows x ncols
    """
    return np.empty((nrows, ncols), dtype=float)


def print_data(filename, a):
    """
    Write data to file in space-delimited format
    :param filename: Name of file to write data to
    :param a: Data that is written to file
    :return: None
    """
    np.savetxt(filename, a, fmt="%16.15f", delimiter=" ")


def calculate_mass(phi, h2, nx, ny):
    ave_mass = np.sum(phi) / (h2 * nx * ny)
    return ave_mass


def f(phi):
    fphi = (1 / 4) * ((phi ^ 2) - 1) ^ 2
    return fphi


def calculate_discrete_energy(phi, h2, nx, ny, epsilon2):
    a = h2 * np.sum(f(phi))
    sum_i = np.sum((phi[2:-1, :] - phi[1:-2, :]) ^ 2)

    b = (epsilon2 / 2) * sum_i
    sum_j = np.sum((phi[:, 2:-1] - phi[:, 1:-2]) ^ 2)

    c = (epsilon2 / 2) * sum_j
    E = a + b + c
    return E


def calculate_discrete_norm_energy(phi, phi0, h2, nx, ny, epsilon2):
    E0 = calculate_discrete_energy(phi0, h2, nx, ny, epsilon2)
    E = calculate_discrete_energy(phi, h2, nx, ny, epsilon2)
    return E / E0
