using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random
using FFTW

include("aux_functions_SAV.jl")
#This function uses the sav method to solve the 
#Cahn-Hilliard equation for the next time step.

#INPUTS
# phi_old  = Prior chemical state.
# r_old    = Prior SAV parameter.
# hx, hy   = Grid spacings in x and y directions.
# k2, k4   = Laplacian and biharmonic operators.
# dt       = Time step.
# epsilon2 = Square of the interface width parameter.
# boundary = 'periodic' or 'neumann' (not explicitly used in this snippet).
# C0       = Regularization parameter.
# Beta     = Relaxation parameter.
# gamma0   = Stabilization parameter.

#OUTPUT
#phi_new = Next chemical state.
#r_new   = Next sav state.
function sav_solver(phi_old, r_old, hx, hy, k2, k4, dt, epsilon2, boundary, C0, Beta, gamma0)

    phi0 = phi_old
    r0 = r_old

    phi0_df = df(phi0, gamma0) #df at phi0
    Lap_dfphi0 = Lap_SAV(phi0_df, k2)    #Lap of df(phi0)
    phi_bar = A_inv_CN(phi0 + dt / 2 * Lap_dfphi0, dt, k2, k4, gamma0, epsilon2)

    # Step 1
    b = b_fun(phi_bar, hx, hy, C0, gamma0)
    g = g_fun_CN(phi0, r0, b, dt, hx, hy, epsilon2, gamma0, Beta, C0, k2)

    AiLb = A_inv_CN(Lap_SAV(b, k2), dt, k2, k4, gamma0, epsilon2)
    Aig = A_inv_CN(g, dt, k2, k4, gamma0, epsilon2)

    gamma = -real.(fft(b .* AiLb))
    gamma = gamma[1, 1] * hx * hy

    # Step 2
    bphi = real.(fft(b .* Aig))
    bphi = bphi[1, 1] * hx * hy / (1 + dt / 4 * gamma)

    # Step 3
    phi_new = dt / 4 * bphi .* AiLb + Aig
    r_new = r_fun(phi_new, phi_old, r0, b, hx, hy, C0, Beta, dt)
    return phi_new, r_new
end