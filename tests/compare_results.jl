#%%
using CSV, DelimitedFiles
using Plots
# or use Makie for interactivity

# --- Step 1: Read in the matrices ---
boundary = "periodic"
timepoint = 200
nx = 128
ny = 128
language1 = "MATLAB"
file1 = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_$(language1)-$(boundary)/$(boundary)_NMG_finaltest"
full_data = readdlm("$(file1)_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

A = full_data[start_row:end_row, 1:nx]

language2 = "MATLAB"
file2 = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_$(language2)-$(boundary)/$(boundary)_SAV_finaltest"
full_data = readdlm("$(file2)_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

B = full_data[start_row:end_row, 1:nx]

# --- Step 2: Compute difference ---
D = A - B

# --- Step 3: Plot heatmaps ---
plot(layout=(2, 3), size=(1200, 800), plot_title="Phi at Timepoint = $(timepoint) (Final)", titlefontsize=16, plot_margin=2)

# Comparing phi
title1 = splitext(basename(file1))[1]
heatmap!(subplot=1, A, title="$(language1) \n $(title1)", colorbar=false, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1), titlefont=font(10),)

title2 = splitext(basename(file2))[1]
heatmap!(subplot=2, B, title="$(language2) \n $(title2)", colorbar=false, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1), titlefont=font(10),)
heatmap!(subplot=3, D, title="Difference: \n $(title1) - $(title2)", colorbar=true, c=:viridis,
    aspect_ratio=1,
    clim=(-maximum(abs.(D)), maximum(abs.(D))), titlefont=font(10),)

# Comparing mass
mass1 = readdlm("$(file1)_mass.csv", ',')
mass2 = readdlm("$(file2)_mass.csv", ',')

mass1 = mass1 .- mass1[1]
mass2 = mass2 .- mass2[1]
plot!(subplot=4, mass1, label="$(language1) \n $(title1)", xlabel="Time Step", ylabel="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))
plot!(subplot=4, mass2, label="$(language2) \n $(title2)", xlabel="Time Step", ylabel="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))

# Comparing energy
energy1 = readdlm("$(file1)_energy.csv", ',')
energy2 = readdlm("$(file2)_energy.csv", ',')
plot!(subplot=5, energy1, label="$(language1) \n $(title1)", xlabel="Time Step", ylabel="Energy", titlefont=font(10), legend=:topright)
plot!(subplot=5, energy2, label="$(language2) \n $(title2)", xlabel="Time Step", ylabel="Energy", titlefont=font(10), legend=:topright)

# L2 norm over time

# calculare L2 norm for each timepoint in phi
l2_norm_err = sqrt.(sum((D) .^ 2, dims=2))
plot!(subplot=6, l2_norm_err, label="$(language1) \n $(title1)", xlabel="Time Step", ylabel="L2 Norm Error", titlefont=font(10), legend=:topright)





