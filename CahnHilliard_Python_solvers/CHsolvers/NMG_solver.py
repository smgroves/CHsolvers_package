import numpy as np
import aux_functions as aux
from laplacian import laplace
from error2 import error2


# note: the C code updates ct, but it doesn't appear to be used anywhere. We'll need to check if it should be returned.
def source(c_old, nx=nx, ny=ny, dt=dt):
    """
    Compute the source term for phi and mu
    :param c_old: phi at a time step
    :return: src_c, the source term for phi, and src_mu, the source term for mu
    """

    src_mu = aux.dmatrix(nx, ny)
    src_c = aux.dmatrix(nx, ny)
    ct = laplace(c_old, nx, ny)
    for i in range(nx):
        for j in range(ny):
            src_c[i, j] = c_old[i, j] / dt - ct[i, j]  # update source term of phi
            src_mu[i, j] = 0  # set source term for mu to zero
    return src_c, src_mu


def relax(
    c_new,
    mu_new,
    su,
    sw,
    nxt,
    nyt,
    c_relax=c_relax,
    xright=xright,
    xleft=xleft,
    dt=dt,
    Cahn=Cahn,
):
    """
    SMOOTH Relaxation operator. This is just solving x =b*A-1 for the system of equations c_new and mu_new, where A is
    the LHS of equations 22 and 23, and b is the RHS.
    :param c_new: c to be smoothed
    :param mu_new: mu to be smoothed
    :param su: sc, locally defined
    :param sw: smu, locally defined
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :param c_relax: number of relaxation operations
    :param xright: right x-coordinate
    :param xleft: left x-coordinate
    :param dt: time step
    :param Cahn: ϵ^2
    :return: c_new, mu_new
    """
    ht2 = ((xright - xleft) / nxt) ** 2  # h2 temp, defined locally
    a = np.empty(4)
    f = np.empty(2)
    # print("c_new before relaxation: \n", c_new)
    # print("mu_new before relaxation: \n", mu_new)
    # print("su before relaxation: \n", su)
    # print("sw before relaxation: \n", sw)
    for iter in range(c_relax):  # c_relax is defined to be 2 in CHsolver.c
        for i in range(nxt):
            for j in range(nyt):
                if i > 0 and i < nxt - 1:
                    x_fac = 2.0
                else:
                    x_fac = 1.0
                if j > 0 and j < nyt - 1:
                    y_fac = 2.0
                else:
                    y_fac = 1.0
                a[0] = 1 / dt
                a[1] = (x_fac + y_fac) / ht2
                a[2] = -(x_fac + y_fac) * Cahn / ht2 - 3 * (c_new[i][j]) ** 2
                a[3] = 1.0

                f[0] = su[i][j]
                f[1] = sw[i][j] - 2 * (
                    c_new[i][j] ** 3
                )  # replaced from c code with a more condensed version
                if (
                    i > 0
                ):  # boundary cases are slightly different because i-1 doesn't exist for i = 0, for example (same for j)
                    f[0] += mu_new[i - 1][j] / ht2
                    f[1] -= Cahn * c_new[i - 1][j] / ht2
                if i < nxt - 1:
                    f[0] += mu_new[i + 1][j] / ht2
                    f[1] -= Cahn * c_new[i + 1][j] / ht2
                if j > 0:
                    f[0] += mu_new[i][j - 1] / ht2
                    f[1] -= Cahn * c_new[i][j - 1] / ht2
                if j < nyt - 1:
                    f[0] += mu_new[i][j + 1] / ht2
                    f[1] -= Cahn * c_new[i][j + 1] / ht2
                det = a[0] * a[3] - a[1] * a[2]
                c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det
                mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det
        # print("f: \n", f)
        # print("a: \n", a)
        # print("c_new: \n", c_new)
        # print("mu_new: \n", mu_new)
    return c_new, mu_new


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


def nonL(c_new, mu_new, nxt, nyt, dt=dt, Cahn=Cahn):
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
            rw[i][j] = mu_new[i][j] - (c_new[i][j]) ** 3 + Cahn * lap_c[i][j]
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


def vcycle(uf_new, wf_new, su, sw, nxf, nyf, ilevel):
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
        Cahn=Cahn,
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
            uf_new=uc_def,
            wf_new=wc_def,
            su=duc,
            sw=dwc,
            nxf=nxc,
            nyf=nyc,
            ilevel=ilevel + 1,
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
            Cahn=Cahn,
        )

    return uf_new, wf_new


def cahn(c_old, c_new, mu, nx=nx, ny=ny, dt=dt, max_it_CH=10000, tol=1e-10):
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx=nx, ny=ny, dt=dt)

    while it_mg2 < max_it_CH and resid2 > tol:

        c_new, mu = vcycle(
            uf_new=c_new, wf_new=mu, su=sc, sw=smu, nxf=nx, nyf=ny, ilevel=1
        )  # TODO why does this give any error when assigned to c_new, mu?
        resid2 = error2(c_old=c_old, c_new=c_new, mu=mu, nxt=nx, nyt=ny, dt=dt)
        it_mg2 += 1

    return c_new
