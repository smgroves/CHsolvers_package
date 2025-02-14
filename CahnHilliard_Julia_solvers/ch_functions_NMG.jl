using LinearAlgebra
using StaticArrays
using Printf

# laplacian function: laplacian(m, nx, ny, h2)
function laplace!(lap_a, a, nxt, nyt, xright, xleft, yright, yleft, boundary)

    ht2 = ((xright - xleft) / nxt)^2
    for i in 1:nxt
        for j in 1:nyt
            if i > 1
                dadx_L = (a[i, j] - a[i-1, j])
            else
                if boundary == "neumann"
                    dadx_L = 0
                elseif boundary == "periodic"
                    dadx_L = a[i, j] - a[nxt-1, j]
                end
            end
            if i < nxt
                dadx_R = (a[i+1, j] - a[i, j])
            else
                if boundary == "neumann"
                    dadx_R = 0
                elseif boundary == "periodic"
                    dadx_R = a[2, j] - a[i, j]
                end
            end
            if j > 1
                dady_B = (a[i, j] - a[i, j-1])
            else
                if boundary == "neumann"
                    dady_B = 0
                elseif boundary == "periodic"
                    dady_B = a[i, j] - a[i, nyt-1]
                end
            end
            if j < nyt
                dady_T = (a[i, j+1] - a[i, j])
            else
                if boundary == "neumann"
                    dady_T = 0
                elseif boundary == "periodic"
                    dady_T = a[i, 2] - a[i, j]
                end
            end
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2
        end
    end
    return lap_a
end

function print_mat(file, matrix)
    open(file, "a", lock=false) do f
        for i = axes(matrix, 1)
            for j = axes(matrix, 2)
                Printf.@printf(f, "%16.15f ", matrix[i, j])
            end
            println(f)
        end
    end
end

function source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary)
    src_mu = zeros(Float64, nx, ny)
    src_c = zeros(Float64, nx, ny)
    ct = zeros(Float64, nx, ny)
    laplace!(ct, c_old, nx, ny, xright, xleft, yright, yleft, boundary)

    src_c .= c_old ./ dt .- ct

    return src_c, src_mu
end

function relax!(c_new, mu_new, su, sw, nxt, nyt, c_relax, xright, xleft, yright, yleft, dt, Cahn)
    ht2 = ((xright - xleft) / nxt)^2
    a = MVector{4,Float64}(undef)
    f = MVector{2,Float64}(undef)
    for iter in 1:c_relax
        for i in 1:nxt
            for j in 1:nyt
                if boundary == "neumann"
                    if i > 1 && i < nxt
                        x_fac = 2.0
                    else
                        x_fac = 1.0
                    end
                    if j > 1 && j < nyt
                        y_fac = 2.0
                    else
                        y_fac = 1.0
                    end
                elseif boundary == "periodic"
                    x_fac = 2.0
                    y_fac = 2.0
                end
                a[1] = 1 / dt
                a[2] = (x_fac + y_fac) / ht2
                a[3] = -(x_fac + y_fac) * Cahn / ht2 - 3 * (c_new[i, j])^2
                cnew_val = (c_new[i, j])
                d2f = -3 * (c_new[i, j])^2

                a[4] = 1.0

                f[1] = su[i, j]
                f[2] = sw[i, j] - 2 * (c_new[i, j])^3

                if i > 1
                    f[1] += mu_new[i-1, j] / ht2
                    f[2] -= Cahn * c_new[i-1, j] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[nxt-1, j] / ht2
                    f[2] -= Cahn * c_new[nxt-1, j] / ht2
                end
                if i < nxt
                    f[1] += mu_new[i+1, j] / ht2
                    f[2] -= Cahn * c_new[i+1, j] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[2, j] / ht2
                    f[2] -= Cahn * c_new[2, j] / ht2
                end
                if j > 1
                    f[1] += mu_new[i, j-1] / ht2
                    f[2] -= Cahn * c_new[i, j-1] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[i, nyt-1] / ht2
                    f[2] -= Cahn * c_new[i, nyt-1] / ht2
                end
                if j < nyt
                    f[1] += mu_new[i, j+1] / ht2
                    f[2] -= Cahn * c_new[i, j+1] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[i, 2] / ht2
                    f[2] -= Cahn * c_new[i, 2] / ht2
                end
                det = a[1] * a[4] - a[2] * a[3]
                c_new[i, j] = (a[4] * f[1] - a[2] * f[2]) / det
                mu_new[i, j] = (-a[3] * f[1] + a[1] * f[2]) / det

            end
        end
    end
    # return c_new, mu_new
end


function restrict_ch(uf, vf, nxc, nyc)
    uc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    vc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    for i in 1:nxc
        for j in 1:nyc
            uc[i, j] = 0.25 * (uf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i - 1), round(Int, 2 * j)] + uf[round(Int, 2 * i), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i), round(Int, 2 * j)])
            vc[i, j] = 0.25 * (vf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i - 1), round(Int, 2 * j)] + vf[round(Int, 2 * i), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i), round(Int, 2 * j)])
        end
    end
    return uc, vc
end

function nonL!(lap_c, lap_mu, c_new, mu_new, nxt, nyt, dt, Cahn, xright, xleft, yright, yleft, boundary)
    ru = zeros(Float64, nxt, nyt)
    rw = zeros(Float64, nxt, nyt)
    laplace!(lap_c, c_new, nxt, nyt, xright, xleft, yright, yleft, boundary)
    laplace!(lap_mu, mu_new, nxt, nyt, xright, xleft, yright, yleft, boundary)

    for i in 1:nxt
        for j in 1:nyt
            ru[i, j] = c_new[i, j] / dt - lap_mu[i, j]
            rw[i, j] = mu_new[i, j] / dt - c_new[i, j]^3 + Cahn * lap_c[i, j]
        end
    end
    return ru, rw
end

# df(c) function: c.^3
# d2f(c) function: 3*c.^2
function defect!(lap_c, lap_mu, uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, Cahn, xright, xleft, yright, yleft)
    ruc, rwc = nonL!(lap_c, lap_mu, uc_new, wc_new, nxc, nyc, dt, Cahn, xright, xleft, yright, yleft, boundary)
    ruf, rwf = nonL!(lap_c, lap_mu, uf_new, wf_new, nxf, nyf, dt, Cahn, xright, xleft, yright, yleft, boundary)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc
end

function prolong_ch(uc, vc, nxc, nyc)
    uf = zeros(Float64, 2 * nxc, 2 * nyc)
    vf = zeros(Float64, 2 * nxc, 2 * nyc)
    for i in 1:nxc
        for j in 1:nyc
            uf[2*i-1, 2*j-1] = uc[i, j]
            uf[2*i-1, 2*j] = uc[i, j]
            uf[2*i, 2*j-1] = uc[i, j]
            uf[2*i, 2*j] = uc[i, j]
            vf[2*i-1, 2*j-1] = vc[i, j]
            vf[2*i-1, 2*j] = vc[i, j]
            vf[2*i, 2*j-1] = vc[i, j]
            vf[2*i, 2*j] = vc[i, j]
        end
    end
    return uf, vf
end

function vcycle(uf_new, wf_new, su, sw, nxf, nyf, ilevel, c_relax, xright, xleft, yright, yleft, dt, Cahn, n_level)
    relax!(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright, xleft, yright, yleft, dt, Cahn)

    if ilevel < n_level
        nxc = trunc(Int64, nxf / 2)
        nyc = trunc(Int64, nyf / 2)
        uc_new, wc_new = restrict_ch(uf_new, wf_new, nxc, nyc)
        lap_c = zeros(Float64, nxf, nyf)
        lap_mu = zeros(Float64, nxf, nyf)
        duc, dwc = defect!(lap_c, lap_mu, uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, Cahn, xright, xleft, yright, yleft)

        uc_def = copy(uc_new)
        wc_def = copy(wc_new)

        uc_def, wc_def = vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1, c_relax, xright, xleft, yright, yleft, dt, Cahn, n_level)


        uc_def = uc_def - uc_new
        wc_def = wc_def - wc_new

        uf_def, wf_def = prolong_ch(uc_def, wc_def, nxc, nyc)

        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def

        relax!(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright, xleft, yright, yleft, dt, Cahn)

    end
    return uf_new, wf_new
end

function error2!(rr, c_old, c_new, mu, nxt, nyt, dt, xright, xleft, yright, yleft, boundary)
    x = 0.0

    rr .= mu .- c_old

    sor = zeros(Float64, nxt, nyt)
    laplace!(sor, rr, nxt, nyt, xright, xleft, yright, yleft, boundary)

    rr .= sor .- (c_new .- c_old) ./ dt


    @simd for i in 1:nxt
        for j in 1:nyt
            x += rr[i, j]^2
        end
    end
    res2 = sqrt(x / (nxt * nyt))
    return res2
end


function cahn!(rr, c_old, c_new, mu, nx, ny, dt, max_it_CH, tol, c_relax, xright, xleft, yright, yleft, Cahn, n_level, boundary; suffix="", print_r=true)
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary)

    while resid2 > tol && it_mg2 < max_it_CH

        c_new, mu = vcycle(c_new, mu, sc, smu, nx, ny, 1, c_relax, xright, xleft, yright, yleft, dt, Cahn, n_level)

        resid2 = error2!(rr, c_old, c_new, mu, nx, ny, dt, xright, xleft, yright, yleft, boundary)
        if print_r
            print_mat("$(outdir)/residual_$(nx)_$(tol)_$(suffix).txt", resid2)
        end
        it_mg2 += 1
    end

    return c_new
end