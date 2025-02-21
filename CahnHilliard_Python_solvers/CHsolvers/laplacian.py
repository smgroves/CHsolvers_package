import aux_functions as aux


def laplace(a, nxt, nyt, xright=xright, xleft=xleft):
    """
    Compute the discrete Laplacian of a
    :param xright:
    :param a: matrix
    :param nxt: nx temp (number of grid points in x-direction, locally defined)
    :param nyt: ny temp (number of grid points in y-direction, locally defined)
    :return: lap_a, the discrete laplacian of a
    """
    lap_a = aux.dmatrix(nxt, nyt)
    h2 = ((xright - xleft) / nxt) ** 2
    for i in range(nxt):
        for j in range(nyt):
            if i > 0:
                dadx_L = a[i, j] - a[i - 1, j]
            else:
                dadx_L = 0
            if i < nxt - 1:
                dadx_R = a[i + 1, j] - a[i, j]
            else:
                dadx_R = 0
            if j > 0:
                dady_B = a[i, j] - a[i, j - 1]
            else:
                dady_B = 0
            if j < nyt - 1:
                dady_T = a[i, j + 1] - a[i, j]
            else:
                dady_T = 0
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / h2
    return lap_a
