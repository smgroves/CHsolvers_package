#%%
using CSV, DelimitedFiles
using Plots
# or use Makie for interactivity

# --- Step 1: Read in the matrices ---
boundary = "neumann"
timepoint = 2
nx = 128
ny = 128

full_data = readdlm("/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_Julia-$(boundary)/$(boundary)_NMG_testv4_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

A = full_data[start_row:end_row, 1:nx]


full_data = readdlm("/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_MATLAB-$(boundary)/$(boundary)_NMG_testv3_phi.csv", ',')

# Extract rows for the desired timepoint
start_row = (timepoint - 1) * ny + 1
end_row = timepoint * ny

B = full_data[start_row:end_row, 1:nx]

# --- Step 2: Compute difference ---
D = A - B

# --- Step 3: Plot heatmaps ---
plot(layout=(1, 3), size=(1200, 400), titlefontsize=10, plot_margin=2)

heatmap!(subplot=1, A, title="Julia: $(boundary) NMG", colorbar=true, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1)
)
heatmap!(subplot=2, B, title="MATLAB: $(boundary) NMG", colorbar=true, c=cgrad([:blue, :white, :red]),
    aspect_ratio=1,
    clim=(-1, 1)
)
heatmap!(subplot=3, D, title="Julia - MATLAB", colorbar=true, c=:viridis,
    aspect_ratio=1,
    clim=(-maximum(abs.(D)), maximum(abs.(D))),
)


