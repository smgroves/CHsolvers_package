using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random

#will need to ensure matrices are floats not ints
function convert_matrix(matrix)
    if eltype(matrix) == Float64
        return matrix  # Matrix is already of type Float64, no need to convert
    elseif eltype(matrix) == Int
        return convert(Matrix{Float64}, matrix)  # Convert to Matrix{Float64}
    else
        throw(ArgumentError("Unsupported matrix element type"))
    end
end

function calculate_mass(phi, h2, nx, ny)
    ave_mass = sum(phi) / (h2 * nx * ny)
    return ave_mass
end

function f(phi)
    fphi = (1 / 4) .* ((phi .^ 2) .- 1) .^ 2
    return fphi
end

function calculate_discrete_energy(phi, h2, nx, ny, Cahn)
    a = h2 * sum(f(phi))
    sum_i = sum((phi[2:end, :] .- phi[1:end-1, :]) .^ 2)

    b = (Cahn / 2) * sum_i
    sum_j = sum((phi[:, 2:end] .- phi[:, 1:end-1]) .^ 2)

    c = (Cahn / 2) * sum_j
    E = a + b + c
    return E
end

function calculate_discrete_norm_energy(phi, phi0, h2, nx, ny, Cahn)
    E0 = calculate_discrete_energy(phi0, h2, nx, ny, Cahn)
    E = calculate_discrete_energy(phi, h2, nx, ny, Cahn)
    return E / E0
end
