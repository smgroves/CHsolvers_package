from . import aux_functions as aux


def initialization(nx, ny):
    """
    Example initialization scheme for phi; this initializes phi for spinodal decomposition
    :return: phi
    """
    phi = aux.dmatrix(nx, ny)
    CPC_width = 5
    cohesin_width = 1
    # for i in range(nx):
    #     for j in range(ny):
    #         phi[i][j] = 0.1 * (1 - 2 * random.randint(0,1))
    for i in range(nx):
        for j in range(ny):
            if i > int(nx / 2 - CPC_width - 1) and i < int(nx / 2 + CPC_width - 1):
                if j > int(ny / 2 - CPC_width - 1) and j < int(ny / 2 + CPC_width - 1):
                    phi[i][j] = 1
                elif i > int(nx / 2 - cohesin_width - 1) and i < int(
                    nx / 2 + cohesin_width - 1
                ):
                    phi[i][j] = 1
                else:
                    phi[i][j] = -1
            else:
                phi[i][j] = -1
    return phi
