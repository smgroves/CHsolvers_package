#%%
using CSV, DelimitedFiles
using Plots
using Printf
# import Pkg;
# Pkg.add("PyPlot");
using PyPlot
# or use Makie for interactivity
pyplot()
# --- Step 1: Read in the matrices ---
boundary = "neumann"
timepoint = 200
nx = 128
ny = 128
language1 = "Julia"
method1 = "NMG"
file1 = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_$(language1)-$(boundary)/$(boundary)_$(method1)_finaltest"
full_data = readdlm("$(file1)_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

A = full_data[start_row:end_row, 1:nx]

language2 = "Julia"
method2 = "SAV"
file2 = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_$(language2)-$(boundary)/$(boundary)_$(method2)_finaltest"
full_data = readdlm("$(file2)_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

B = full_data[start_row:end_row, 1:nx]

# --- Step 2: Compute difference ---
D = A - B

# --- Step 3: Plot heatmaps ---
plot(layout=(2, 3), size=(1200, 800), plot_title="Phi at Timepoint = $(timepoint) (Final)", titlefontsize=16)

# Comparing phi
title1 = splitext(basename(file1))[1]
heatmap!(subplot=1, A, title="$(language1) \n $(title1)", colorbar=false, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1), titlefont=font(10), interpolate=false,      # key to prevent smoothing/bleeding
    framestyle=:box,)

title2 = splitext(basename(file2))[1]
heatmap!(subplot=2, B, title="$(language2) \n $(title2)", colorbar=false, c=cgrad([:blue, :white, :red]),
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
plot!(subplot=4, mass1, label="$(language1) \n $(title1)", xlabel="Time Step", title="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))
plot!(subplot=4, mass2, label="$(language2) \n $(title2)", xlabel="Time Step", title="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))

# Comparing energy
energy1 = readdlm("$(file1)_energy.csv", ',')
energy2 = readdlm("$(file2)_energy.csv", ',')
plot!(subplot=5, energy1, label="$(language1) \n $(title1)", xlabel="Time Step", title="Energy", titlefont=font(10), legend=:topright)
plot!(subplot=5, energy2, label="$(language2) \n $(title2)", xlabel="Time Step", title="Energy", titlefont=font(10), legend=:topright)

# L2 norm over time

# calculare L2 norm for each timepoint in phi
l2_norm_err = sqrt.(sum((D) .^ 2, dims=2))
plot!(subplot=6, l2_norm_err, label="Error", xlabel="Time Step", title="L2 Norm Error", titlefont=font(10), legend=:topright)





