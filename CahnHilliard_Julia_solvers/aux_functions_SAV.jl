using FFTW

function ext(x) #good
    nx, ny = size(x)
    x_ext = zeros(2 * nx, 2 * ny)

    # Original block
    x_ext[1:nx, 1:ny] = x

    # Flip horizontally
    x_ext[1:nx, ny+1:2*ny] = x[:, end:-1:1]

    # Flip vertically
    x_ext[nx+1:2*nx, 1:ny] = x[end:-1:1, :]

    # Flip both
    x_ext[nx+1:2*nx, ny+1:2*ny] = x[end:-1:1, end:-1:1]

    return x_ext
end

function f(phi) #good
    # f = @(x) 0.25*(x.^2-1).^2
    fphi = (phi .^ 2 .- 1) .^ 2 ./ 4
    return fphi
end

function r0_fun(phi0, hx, hy, C0) #good
    fphi = f(phi0)
    r0 = sqrt(hx * hy * sum(fphi) + C0)
    return r0
end

function extback(x_ext) #good
    # Shrinks from 2*nx x 2*ny back to nx x ny (upper-left block)
    nx_ext, ny_ext = size(x_ext)
    nx = nx_ext ÷ 2
    ny = ny_ext ÷ 2
    x_back = x_ext[1:nx, 1:ny]

    return x_back
end

function df(phi, gamma0) #good
    return phi .^ 3 .- (1 .+ gamma0) .* phi
end

function Lap_SAV(phi, k2) #good
    return real.(ifft(k2 .* fft_filtered(phi)))
end

function A_inv_CN(phi, dt, k2, k4, gamma0, epsilon2) #good
    denom = 1 .+ (dt / 2) * epsilon2 .* k4 .- (dt / 2) * gamma0 .* k2
    return real.(ifft(fft_filtered(phi) ./ denom))
end

function fft_filtered(x) #removed real(), good
    return (fft(x))
end

function b_fun(phi, hx, hy, C0, gamma0) #good
    E1 = fft_filtered(f(phi))
    return df(phi, gamma0) ./ sqrt.(E1[1, 1] * hx * hy .+ C0)
end

function g_fun_CN(phi0, r0, b, dt, hx, hy, epsilon2, gamma0, Beta, C0, k2) #good
    Lap_phi0 = Lap_SAV(phi0, k2)
    Lap_Lap_phi0 = Lap_SAV(Lap_phi0, k2)

    bphi0 = fft_filtered(b .* phi0)
    bphi0 = hx * hy * bphi0[1, 1]

    E1 = fft_filtered(f(phi0))

    g = phi0 .- (dt / 2) * epsilon2 .* Lap_Lap_phi0 .+ (dt / 2) * gamma0 .* Lap_phi0 .+
        dt .* Lap_SAV(b, k2) .* (r0 .- (1 / 4) * bphi0 .- (1 / 2) * Beta * dt * r0 .* (r0 .- sqrt.(E1[1, 1] * hx * hy .+ C0)))

    return g
end

function r_fun(phi, phi0, r0, b, hx, hy, C0, Beta, dt)
    bphi0 = fft_filtered(b .* phi0)
    bphi0 = hx * hy * bphi0[1, 1]

    bphi = fft_filtered(b .* phi)
    bphi = hx * hy * bphi[1, 1]

    E1 = fft_filtered(f(phi0))

    r = r0 .+ (1 / 2) * (bphi .- bphi0) .- Beta * dt * r0 .* (r0 .- sqrt.(E1[1, 1] * hx * hy .+ C0))

    return r
end

