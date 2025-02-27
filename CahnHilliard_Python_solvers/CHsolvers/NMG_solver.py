import numpy as np
import aux_functions as aux
from laplacian import laplace
from error2 import error2


# note: the C code updates ct, but it doesn't appear to be used anywhere. We'll need to check if it should be returned.
def source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary):
    """
    Compute the source term for phi and mu
    :param c_old: phi at a time step
    :return: src_c, the source term for phi, and src_mu, the source term for mu
    """

    src_mu = aux.dmatrix(nx, ny)
    src_c = aux.dmatrix(nx, ny)
    ct = laplace(c_old, nx, ny, xright, xleft, yright, yleft, boundary)
    for i in range(nx):
        for j in range(ny):
            src_c[i, j] = c_old[i, j] / dt - ct[i, j]  # update source term of phi
            src_mu[i, j] = 0  # set source term for mu to zero
    return src_c, src_mu


def restrict_ch(uf, vf, nxc, nyc):
    """
    Restrict the defect twofold in each direction
    uf and vf get compressed to uc and vc with dimensions nxc and nyc
    Note that changing from C to Python requires adding 1 instead of subtracting in formulas
    :param uf: uf matrix to be restricted
    :param vf: vf matrix to be restricted
    :param nxc: number of grid points in x-direction of uc
    :param nyc: number of grid points in y-direction of vc
    :return: uc, vc
    """
    uc = aux.dmatrix(nxc, nyc)
    vc = aux.dmatrix(nxc, nyc)
    for i in range(nxc):
        for j in range(nyc):
            uc[i][j] = 0.25 * (
                uf[2 * i][2 * j]
                + uf[2 * i + 1][2 * j]
                + uf[2 * i][2 * j + 1]
                + uf[2 * i + 1][2 * j + 1]
            )
            vc[i][j] = 0.25 * (
                vf[2 * i][2 * j]
                + vf[2 * i + 1][2 * j]
                + vf[2 * i][2 * j + 1]
                + vf[2 * i + 1][2 * j + 1]
            )
    return uc, vc


def nonL(c_new, mu_new, nxt, nyt, dt, epsilon2):
    """
    NSO operator
    :param c_new: c at a time step
    :param mu_new: mu at a time step
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :return: ru, rw
    """
    ru = aux.dmatrix(nxt, nyt)
    rw = aux.dmatrix(nxt, nyt)

    lap_c = laplace(c_new, nxt, nyt)
    lap_mu = laplace(mu_new, nxt, nyt)
    for i in range(nxt):
        for j in range(nyt):
            ru[i][j] = c_new[i][j] / dt - lap_mu[i][j]
            rw[i][j] = mu_new[i][j] - (c_new[i][j]) ** 3 + epsilon2 * lap_c[i][j]
    return ru, rw


def defect(uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc):
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc


def prolong_ch(uc, vc, nxc, nyc):
    uf = np.zeros((2 * nxc, 2 * nyc))
    vf = np.zeros((2 * nxc, 2 * nyc))
    for i in range(nxc):
        for j in range(nyc):
            uf[2 * i][2 * j] = uf[2 * i + 1][2 * j] = uf[2 * i][2 * j + 1] = uf[
                2 * i + 1
            ][2 * j + 1] = uc[i][j]
            vf[2 * i][2 * j] = vf[2 * i + 1][2 * j] = vf[2 * i][2 * j + 1] = vf[
                2 * i + 1
            ][2 * j + 1] = vc[i][j]
    return uf, vf


def vcycle(
    uf_new,
    wf_new,
    su,
    sw,
    nxf,
    nyf,
    ilevel,
    c_relax,
    xright,
    xleft,
    yright,
    yleft,
    dt,
    epsilon2,
    n_level,
    boundary,
):
    """
    FAS multigrid cycle
    """
    # relax the input data
    # print("before relaxing IN VCYCLE")
    # print("nxf: ", nxf) #same
    # print("uf_new: \n", uf_new) #same
    # print("wf_new: \n", wf_new)
    # print("su: \n", su)
    # print("sw: \n", sw)
    uf_new, wf_new = relax(
        c_new=uf_new,
        mu_new=wf_new,
        su=su,
        sw=sw,
        nxt=nxf,
        nyt=nyf,
        c_relax=c_relax,
        xright=xright,
        xleft=xleft,
        dt=dt,
        epsilon2=epsilon2,
    )
    # print("after relaxing IN VCYCLE")
    # print("uf_new: \n", uf_new) #same
    # print("wf_new: \n", wf_new)
    # print("su: \n", su)
    # print("sw: \n", sw)
    # If the number of multigrid levels has not been reached
    if ilevel < n_level:
        # print("ilevel: ", ilevel)
        nxc = int(nxf / 2)
        nyc = int(nyf / 2)
        uc_new, wc_new = restrict_ch(uf=uf_new, vf=wf_new, nxc=nxc, nyc=nyc)

        duc, dwc = defect(
            uf_new=uf_new,
            wf_new=wf_new,
            suf=su,
            swf=sw,
            nxf=nxf,
            nyf=nyf,
            uc_new=uc_new,
            wc_new=wc_new,
            nxc=nxc,
            nyc=nyc,
        )

        uc_def = uc_new.copy()
        wc_def = wc_new.copy()

        uc_def, wc_def = vcycle(
            uc_def,
            wc_def,
            duc,
            dwc,
            nxc,
            nyc,
            ilevel + 1,
            c_relax,
            xright,
            xleft,
            dt,
            epsilon2,
            n_level,
            boundary,
        )

        uc_def = uc_def - uc_new
        wc_def = wc_def - wc_new

        uf_def, wf_def = prolong_ch(uc=uc_def, vc=wc_def, nxc=nxc, nyc=nyc)

        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def

        uf_new, wf_new = relax(
            c_new=uf_new,
            mu_new=wf_new,
            su=su,
            sw=sw,
            nxt=nxf,
            nyt=nyf,
            c_relax=c_relax,
            xright=xright,
            xleft=xleft,
            dt=dt,
            epsilon2=epsilon2,
        )

    return uf_new, wf_new


def cahn(c_old, c_new, mu, nx, ny, dt, solver_iter, tol, xright, xleft):
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx, ny, dt, xright, xleft)

    while it_mg2 < solver_iter and resid2 > tol:

        c_new, mu = vcycle(
            uf_new=c_new, wf_new=mu, su=sc, sw=smu, nxf=nx, nyf=ny, ilevel=1
        )
        resid2 = error2(c_old=c_old, c_new=c_new, mu=mu, nxt=nx, nyt=ny, dt=dt)
        it_mg2 += 1

    return c_new
