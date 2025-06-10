using DelimitedFiles
using Dates
using BenchmarkTools
using Printf
include("CahnHilliard_NMG.jl")
include("CahnHilliard_SAV.jl")
include("ch_movie_from_file.jl")

indir = "./IC/"
outdir = "./output/output_julia/"
n_relax = 4
GridSize = 64
h = 1 / GridSize
m = 8
epsilon = m * h / (2 * sqrt(2) * atanh(0.9))

boundary = "periodic"
printphi = true
dt_out = 1
dt = 5.5e-6
max_it = 10

phi0 = zeros(GridSize, GridSize)
phi0[1, 10] = 2
phi0[32, 1] = 6
phi0[50, 64] = 8
phi0[32, 32] = 10
phi0[64, 25] = 4

################
# NMG
################
method = "NMG"
pathname = @sprintf("%speriodic_NMG_testv2_", outdir)
result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_NMG(phi0; t_iter=max_it, dt=dt, dt_out=dt_out, m=m, boundary=boundary, printphi=printphi, pathname=pathname)

t_out = result[1]
phi_t = result[2]
mass_t = result[3]
E_t = result[4]

open("$(pathname)mass.csv", "w", lock=false) do f
    writedlm(f, mass_t, " ")
end
open("$(pathname)energy.csv", "w", lock=false) do f
    writedlm(f, E_t, " ")
end

########
# SAV
########
# method = "SAV"
# pathname = @sprintf("%speriodic_SAV_test_", outdir)
# result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_SAV(phi0; t_iter=max_it, dt=dt, dt_out=dt_out, m=m, boundary=boundary, printphi=printphi, pathname=pathname)


# t_out = result[1]
# phi_t = result[2]
# mass_t = result[3]
# E_t = result[4]

# open("$(pathname)mass.csv", "w", lock=false) do f
#     writedlm(f, mass_t, " ")
# end
# open("$(pathname)energy.csv", "w", lock=false) do f
#     writedlm(f, E_t, " ")
# end
