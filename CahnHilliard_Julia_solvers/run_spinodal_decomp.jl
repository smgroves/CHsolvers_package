using DelimitedFiles
using Dates
using BenchmarkTools
include("CahnHilliard_NMG.jl")
include("ch_movie_from_file.jl")
indir = "./IC/"
outdir = "./output/output_julia/"
################################
# Spinodal decomposition smoothed NMG
################################
method = "NMG"
n_relax = 4
for GridSize in [128]
    h = 1 / GridSize
    m = 8
    epsilon = m * h / (2 * sqrt(2) * atanh(0.9))

    dt = 1E-05;
    max_it = 20000;
    boundary = "neumann"
    printphi = true
    phi0 = initialization_from_file("$(indir)initial_phi_$(GridSize)_smooth_n_relax_$n_relax.csv", GridSize, GridSize, delim=',', transpose_matrix=false)
    pathname = @sprintf("%s%s_Julia_%d_dt_%.2e_Nx_%d_n_relax_%d_", outdir, method, max_it, dt, GridSize, n_relax);
    dt_out = 10
    date_time = now()
    result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_NMG(phi0; t_iter=max_it, dt=dt, dt_out=dt_out, m=m, boundary=boundary, printphi=printphi, pathname=pathname)
    open("./Job_specs.csv", "a", lock=false) do f
        writedlm(f, [Dates.format(date_time, "mm/dd/yyyy HH:MM:SS") "NMG_spinodal_decomp_smoothed_dtout_2000" "Julia" "NMG" GridSize epsilon dt 1e-5 max_it 1e4 pathname elapsed_time mem_allocated / 1e6], ",")
    end

    t_out = result[1]
    phi_t = result[2]
    mass_t = result[3]
    E_t = result[4]
    # open("$(pathname)t_out.csv", "w", lock=false) do f
    #     writedlm(f, t_out, " ")
    # end
    # open("$(pathname)mass.csv", "w", lock=false) do f
    #     writedlm(f, mass_t, " ")
    # end
    # open("$(pathname)energy.csv", "w", lock=false) do f
    #     writedlm(f, E_t, " ")
    # end
    # t_out = readdlm("$(pathname)t_out.csv", ',', Float64)
    ch_movie_from_file("$(pathname)phi.csv", t_out, 128; dtframes=10, filename="$(pathname)movie_jl", filetype="mp4", colorbar_type="default")
end