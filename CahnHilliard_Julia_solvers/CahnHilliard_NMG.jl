# Julia implementation of a Multigrid solver for the Cahn-Hilliard equation with various optimizations for performance (particularly memory allocation)
# Author: Sarah Groves
# October 8, 2024

#=================================================
# Summary of Optimizations:
# [ ] Preallocate arrays and pass them to functions.
# [ ] Use in-place operations and @. for broadcasting.
# [ ] Add @inbounds and @simd for loop optimization.
# [ ] Use StaticArrays for small fixed-size arrays.
# [ ] Replace global constants with function arguments.
# [ ] Use @views to avoid array copying.
# [ ] Profile and benchmark to identify bottlenecks.
# [ ] Access arrays in a column-major order.
# [ ] Precompute and reuse frequently calculated constants.
# [ ] Consider SparseArrays for sparse matrices.
# [ ] Ensure type stability with type annotations.
=================================================#

using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random
include("ch_functions_NMG.jl")
include("ch_initialization.jl")

Random.seed!(1234) #note that when using a random seed, you must RESTART the REPL or run in a new instance for the runs to be the same. 


function multigrid_solver(oc, nx, tol, outdir; max_it=1000, max_it_CH=10000, suffix="", overwrite=true, print_phi=true, print_mass=false, print_e=false, print_r=true, dt=2.5e-5, m=8, ns=10, epsilon=0.0, check_dir=true, c_relax=2, domain=[1.0, 0.0, 1.0, 0.0], boundary="neumann")
    while true
        xright, xleft, yright, yleft = domain
        ny = nx
        n_level::Int = trunc(log(nx) / log(2.0) + 0.1)  # original c code uses natural log too
        h = xright / nx  # space step size defined as a global variable
        h2 = h^2 #space step size squared defined as a global variable

        if epsilon == 0
            epsilon = m * h / (2 * sqrt(2) * atanh(0.9))
        end
        Cahn = epsilon^2  # ϵ^2 defined as a global variable

        #make directory
        if isdir(outdir)
            if !isempty(outdir)
                if overwrite == false
                    if check_dir == true
                        println("Warning: Directory is not empty. Results may be appended to existing files. Are you sure you want to continue? (Y/N)")
                        input = readline()
                        if input == "Y" || input == "y"
                            println("Appending to any existing files.")
                        else
                            println("End.")
                            break
                        end
                    end
                else
                    if check_dir == true
                        println("Warning: overwriting directory with new files. Are you sure you want to continue? (Y/N)")
                        input = readline()
                        if input == "Y" || input == "y"
                            rm(outdir, recursive=true)
                            mkdir(outdir)
                        else
                            println("End.")
                            break
                        end
                    end
                end
            end
        else
            mkpath(outdir)
        end

        #check if oc is in the range -1 to 1; if not, rescale to 0 to 1 and then shift with 2x - 1
        if minimum(oc) >= 0
            type = "psi"
            oc = (oc .- minimum(oc)) ./ (maximum(oc) - minimum(oc))
            oc = 2 .* oc .- 1    # psi0=(phi0+1)/2 
        else #if minimum(oc)>=-1 && maximum(oc)<=1
            type = "phi"
            # else
            # println("Warning: cannot interpret initial data. Data should be in range [-1,1] or only positive values.")
            # break
        end
        println(type)

        #initialization
        println("nx = $nx, ny = $ny, dt = $dt, epsilon = $epsilon, max_it = $max_it,max_it_CH= $max_it_CH, ns = $ns, n_level = $n_level")
        mu = zeros(Float64, nx, ny)
        nc = copy(oc)
        oc0 = copy(oc)
        if print_phi
            if type == "phi"
                open("$(outdir)/phi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "w", lock=false) do f
                    writedlm(f, oc, " ")
                end
            elseif type == "psi"
                open("$(outdir)/psi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "w", lock=false) do f
                    psi = (oc .+ 1) ./ 2
                    writedlm(f, psi, " ")
                end
            end
        end
        rr = zeros(Float64, nx, ny)
        #run cahn hilliard
        for it in 1:max_it
            cahn!(rr, oc, nc, mu, nx, ny, dt, max_it_CH, tol, c_relax, xright, xleft, yright, yleft, Cahn, n_level, boundary, suffix=suffix, print_r=print_r)

            if print_mass
                print_mat("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_mass(oc, h2, nx, ny))
            end
            if print_e
                print_mat("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_discrete_norm_energy(oc, oc0, h2, nx, ny, Cahn))
            end
            oc = copy(nc)
            if it % ns == 0
                if print_phi
                    println(it)
                    if type == "phi"
                        open("$(outdir)/phi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "a", lock=false) do f
                            writedlm(f, oc, " ")
                        end
                    elseif type == "psi"
                        open("$(outdir)/psi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "a", lock=false) do f
                            psi = (oc .+ 1) ./ 2
                            writedlm(f, psi, " ")
                        end
                    end
                end
            end
        end
        break
    end
end

function CahnHilliard_NMG()

end