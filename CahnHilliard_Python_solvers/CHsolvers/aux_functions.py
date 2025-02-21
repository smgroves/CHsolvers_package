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
