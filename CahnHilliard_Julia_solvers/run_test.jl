# # FIGURE 1
using DelimitedFiles
using Dates
using BenchmarkTools
using Printf
include("CahnHilliard_NMG.jl")
include("CahnHilliard_SAV.jl")
include("ch_movie_from_file.jl")
indir = "./IC/"
boundary = "neumann"

outdir = @sprintf("./output/output_Julia-%s", boundary);

m = 8;
GridSize = 128;
h = 1 / GridSize;
epsilon = m * h / (2 * sqrt(2) * atanh(0.9));
dt = 5.5e-6;
max_it = 100;

print_phi = true;
dt_out = 1;
nx = GridSize;
ny = GridSize;
n_relax = 4;
# #################################################
init_file = @sprintf("%s/initial_phi_%d_smooth_n_relax_%d.csv", indir, GridSize, n_relax);
phi0 = initialization_from_file(init_file, GridSize, GridSize, delim=',', transpose_matrix=false)

###################################################
##RUN NMG SOLVER 
###################################################

pathname = @sprintf("%s/%s_NMG_testv3_", outdir, boundary);
println("Running NMG solver with parameters: %s\n", pathname);
result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_NMG(phi0,
    t_iter=max_it,
    dt=dt,
    m=m,
    boundary=boundary,
    printphi=print_phi,
    pathname=pathname,
    dt_out=dt_out);

t_out = result[1]
phi_t = result[2]
mass_t = result[3]
E_t = result[4]
open("$(pathname)t_out.csv", "w", lock=false) do f
    writedlm(f, t_out, " ")
end
open("$(pathname)mass.csv", "w", lock=false) do f
    writedlm(f, mass_t, " ")
end
open("$(pathname)energy.csv", "w", lock=false) do f
    writedlm(f, E_t, " ")
end




###################################################
##RUN SAV SOLVER 
###################################################

pathname = @sprintf("%s/%s_SAV_testv3_", outdir, boundary);
println("Running SAV solver with parameters: %s\n", pathname);
result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_SAV(phi0,
    t_iter=max_it,
    dt=dt,
    m=m,
    boundary=boundary,
    printphi=print_phi,
    pathname=pathname,
    dt_out=dt_out);

t_out = result[1]
phi_t = result[2]
mass_t = result[3]
E_t = result[4]
open("$(pathname)t_out.csv", "w", lock=false) do f
    writedlm(f, t_out, " ")
end
open("$(pathname)mass.csv", "w", lock=false) do f
    writedlm(f, mass_t, " ")
end
open("$(pathname)energy.csv", "w", lock=false) do f
    writedlm(f, E_t, " ")
end
