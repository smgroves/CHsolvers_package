#%%
using CSV, DelimitedFiles
using Plots
using Printf
using Statistics
# import Pkg;
# Pkg.add("PyPlot");
using PyPlot
# or use Makie for interactivity
pyplot()
# --- Step 1: Read in the matrices ---
boundary = "periodic"

folder = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_julia/n_relax_test/$(boundary)"
n_relax = 2
timepoint = 100
nx = 128
ny = 128
method1 = "NMG"
file1 = "$(folder)/$(method1)_$(boundary)_Julia_2000_dt_5.50e-06_Nx_128_n_relax_$(n_relax)"
full_data1 = readdlm("$(file1)_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

A = full_data1[start_row:end_row, 1:nx]

method2 = "SAV"
file2 = "$(folder)/$(method2)_$(boundary)_Julia_2000_dt_5.50e-06_Nx_128_n_relax_$(n_relax)"
full_data2 = readdlm("$(file2)_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny
println(size(full_data2))
B = full_data2[start_row:end_row, 1:nx]

# --- Step 2: Compute difference ---
D = A - B

# --- Step 3: Plot heatmaps ---
println("Plotting heatmaps...")
plot(layout=(2, 3), size=(1200, 800), plot_title="", titlefontsize=16)

# Comparing phi
title1 = splitext(basename(file1))[1]
heatmap!(subplot=1, A, title="$(title1)\nFinal timepoint", colorbar=false, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1), titlefont=font(10), interpolate=false,      # key to prevent smoothing/bleeding
    framestyle=:box,)

title2 = splitext(basename(file2))[1]
heatmap!(subplot=2, B, title="$(title2)\nFinal timepoint", colorbar=false, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1), titlefont=font(10), interpolate=false,      # key to prevent smoothing/bleeding
    framestyle=:box,)

maxval = maximum(abs.(D))
exponent = floor(Int, log10(maxval))
scale = 10.0^exponent
D_scaled = D ./ scale

# Step 2: set ticks around scaled data
# ticks = round.(range(-1, 1; length=5), digits=1)
low = -(maximum(abs.(D_scaled)))
high = -low
ticks = round.(range(low, high; length=5), digits=2)
ticklabels = [Printf.@sprintf("%.1fe%d", t, exponent) for t in ticks]

# Step 3: plot
heatmap!(subplot=3, D_scaled, title="Difference: \n $(title1) - $(title2)",
    clim=(-maximum(abs.(D_scaled)), maximum(abs.(D_scaled))),
    colorbar_ticks=(ticks, ticklabels),
    c=:viridis,
    aspect_ratio=1,
    titlefont=font(10),
    interpolate=false,      # key to prevent smoothing/bleeding
    framestyle=:box,
)
# heatmap!(D, colorbar=true, c=:viridis,
#     aspect_ratio=1, colorbar_tickformatter=x -> Printf.sprintf("%.1e", x),
#     clim=(-maximum(abs.(D)), maximum(abs.(D))), titlefont=font(10), aect_ratio=1
# )

# Comparing mass
mass1 = readdlm("$(file1)_mass.csv", ',')
mass2 = readdlm("$(file2)_mass.csv", ',')

mass1 = mass1 .- mass1[1]
mass2 = mass2 .- mass2[1]
plot!(subplot=4, mass1, label="$(title1)", xlabel="Time Step", title="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))
plot!(subplot=4, mass2, label="$(title2)", xlabel="Time Step", title="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))

# Comparing energy
energy1 = readdlm("$(file1)_energy.csv", ',')
energy2 = readdlm("$(file2)_energy.csv", ',')
plot!(subplot=5, energy1, label="$(title1)", xlabel="Time Step", title="Energy", titlefont=font(10), legend=:topright)
plot!(subplot=5, energy2, label="$(title2)", xlabel="Time Step", title="Energy", titlefont=font(10), legend=:topright)

# L2 norm over time

# calculare L2 norm for each timepoint in phi
l2_norm_err = sqrt.(sum((D) .^ 2, dims=2))
ave_err = mean(l2_norm_err)
plot!(subplot=6, l2_norm_err, label="Error", xlabel="Time Step", title="L2 Norm Error", titlefont=font(10), legend=:topright)
hline!(subplot=6, [ave_err], linestyle=:dot, color=:black, label="Average Error = $(round(ave_err, digits = 4))")  # add horizontal dotted line at y = 0.0


#%%
# for boundary = ["neumann", "periodic"]
ave_err = Vector{Float64}()
push!(ave_err, 0.2)
for n_relax = [2, 4, 8, 16]
    println(n_relax)
    # --- Step 1: Read in the matrices ---
    boundary = "periodic"
    folder = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_julia/n_relax_test/$(boundary)"
    timepoint = 100
    nx = 128
    ny = 128
    method1 = "NMG"
    file1 = "$(folder)/$(method1)_$(boundary)_Julia_2000_dt_5.50e-06_Nx_128_n_relax_$(n_relax)"
    full_data1 = readdlm("$(file1)_phi.csv", ',')

    # Extract rows for the desired timepoint
    start_row = (timepoint - 1) * ny + 1
    end_row = timepoint * ny

    A = full_data1[start_row:end_row, 1:nx]

    method2 = "SAV"
    file2 = "$(folder)/$(method2)_$(boundary)_Julia_2000_dt_5.50e-06_Nx_128_n_relax_$(n_relax)"
    full_data2 = readdlm("$(file2)_phi.csv", ',')

    # Extract rows for the desired timepoint
    start_row = (timepoint - 1) * ny + 1
    end_row = timepoint * ny
    B = full_data2[start_row:end_row, 1:nx]

    # --- Step 2: Compute difference ---
    D = A - B

    l2_norm_err = sqrt.(sum((D) .^ 2, dims=2))
    #append to ave_err
    println(mean(l2_norm_err))
    push!(ave_err, mean(l2_norm_err))

end
print(ave_err)
#plot ave_err vs n_relax
n_relax = [0, 2, 4, 8, 16]
plot(n_relax, ave_err, marker=:o, xlabel="n_relax", ylabel="Average L2 Norm Error", title="Average L2 Norm Error vs n_relax", titlefont=font(10), legend=false)



