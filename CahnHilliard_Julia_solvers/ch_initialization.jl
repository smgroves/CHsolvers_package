using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random

function initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
    phi = zeros(Float64, nx, ny)
    # CPC_width = 20
    # cohesin_width = 4
    @simd for i in 1:nx
        for j in 1:ny
            if i > round(nx / 2) - CPC_width && i < round(nx / 2) + CPC_width
                if j > round(ny / 2) - CPC_width && j < round(ny / 2) + CPC_width
                    phi[i, j] = 1
                elseif i > round(nx / 2) - cohesin_width && i < round(nx / 2) + cohesin_width
                    phi[i, j] = 1
                else
                    phi[i, j] = -1
                end
            else
                phi[i, j] = -1
            end
        end
    end
    return phi
end

function initialize_round_CPC(nx, ny; CPC_width=10, cohesin_width=4)
    # Create an empty matrix filled with -1
    phi = fill(-1.0, nx, ny)

    # Define the center of the matrix
    center = (nx) / 2

    # Loop through each element of the matrix
    @simd for i in 1:nx
        for j in 1:ny
            # Calculate the distance from the center
            distance = norm([i - center, j - center])

            # Check if the distance is less than or equal to CPC_width
            if distance <= CPC_width
                phi[i, j] = 1.0
            elseif i > round((nx) / 2) - cohesin_width && i < round((nx) / 2) + cohesin_width
                phi[i, j] = 1.0
            end
        end
    end
    return phi
end


function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function initialization_from_function(nx, ny, h; R0=0.1, epsilon=0.01)
    phi = zeros(Float64, nx, ny)
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    R = @.sqrt((xx - 0.5)^2 + (yy - 0.5)^2)
    # eps_c = epsilon
    # delta = eps_c * sqrt(2)
    # psi0 = 0.5 * (1 .+ @.tanh((R0 .- R) / (2 * delta)))
    # phi = 2 .* psi0 .- 1    # psi0=(phi0+1)/2
    phi = @.tanh((R0 .- R) / (sqrt(2) * epsilon))
    return phi
end

function initialization_random(nx, ny)
    return (1 .- 2 .* rand(nx, ny))
end

function initialization_spinodal(nx, ny)
    return (rand([-1.0, 1.0], nx, ny))
end

function initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
    phi = readdlm(file, delim, Float64)
    if size(phi) != (nx, ny)
        print("Warning: phi from file is wrong size: $(size(phi)) Expected: $(nx), $(ny)")
    end
    if transpose_matrix
        phi = transpose(phi)
    end
    return phi
end

function ch_initialization(nx, ny; method="spinodal", initial_file="", delim=",", h=1 / 128, R0=0.1, epsilon=0.01, cohesin_width=4, CPC_width=20)
    if method == "random"
        phi0 = initialization_random(nx, ny)
    elseif method == "droplet"
        phi0 = initialization_from_function(nx, ny, h, R0=R0, epsilon=epsilon)
    elseif method == "geometric"
        phi0 = initialize_geometric_CPC(nx, ny, CPC_width=CPC_width, cohesin_width=cohesin_width)
    elseif method == "file"
        phi0 = initialization_from_file(initial_file, nx, ny, delim=delim)
    elseif method == "spinodal"
        phi0 = initialization_spinodal(nx, ny)
    else
        println("Warning: initialize must be one of [random, droplet, geometric, file, spinodal].")
    end
    return phi0
end
