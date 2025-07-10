#%%
println("Starting...")
using CSV, DelimitedFiles
using Plots
using Printf
using Statistics
# import Pkg;
# Pkg.add("PyPlot");
# using PyPlot
# or use Makie for interactivity
# pyplot()
# --- Step 1: Read in the matrices ---
boundary = "neumann"

#%%
folder = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_rivanna_07_2025"
nx = 512
file1 = "NMG_Julia_2000_dt_5.50e-06_Nx_512_neumann_dtout_1050p_phi"
full_data1 = readdlm("$(folder)/$(file1).csv", ',')

# Extract rows for the desired timepoint
# start_row = (timepoint - 1) * nx + 1
# end_row = timepoint * nx

# A = full_data1[start_row:end_row, 1:nx]
#%%
file2 = "SAV_Julia_2000_dt_5.50e-06_Nx_512_neumann_dtout_1050p_phi"
full_data2 = readdlm("$(folder)/$(file2).csv", ',')

#%%

#reshape full_data1 from [AxB,C] to [A,B,C]
timepoints = 201  # number of timepoints
tp1 = size(full_data1, 1) ÷ nx
tp2 = size(full_data2, 1) ÷ nx
if tp1 != tp2
    println("Timepoints do not match: $(tp1) vs $(tp2)")
end
full_data1 = reshape(full_data1, (nx, nx, tp1))
full_data1 = full_data1[:, :, 1:min(timepoints, tp2, tp1)]  # ensure we only take the first min(timepoints, tp1) timepoints
full_data2 = reshape(full_data2, (nx, nx, tp2))
full_data2 = full_data2[:, :, 1:min(timepoints, tp1, tp2)]  # ensure we only take the first min(timepoints, tp2) timepoints
# full_data2 = full_data2[:, :, 2:202]

#%%

# calculare L2 norm for each timepoint in phi
rmse = vec(sqrt.(mean((full_data1 - full_data2) .^ 2, dims=(1, 2))))
ave_err = mean(rmse)
plot(rmse, label="RMSE", xlabel="Time Step", title="RMSE \n$(file1)\n$(file2)", titlefont=font(10), legend=:topleft, xlims=(1, min(timepoints, tp1, tp2)))
hline!([ave_err], linestyle=:dot, color=:black, label="Average RMSE = $(round(ave_err, digits = 4))")  # add horizontal dotted line at y = 0.0

#save FIGURE
savefig("/Users/smgroves/Documents/GitHub/CHsolvers_package/tests/rivanna_output/$(file1)_$(file2).png")

# # Comparing mass
# mass1 = readdlm("$(file1)_mass.csv", ',')
# mass2 = readdlm("$(file2)_mass.csv", ',')

# mass1 = mass1 .- mass1[1]
# mass2 = mass2 .- mass2[1]
# plot!(subplot=4, mass1, label="$(title1)", xlabel="Time Step", title="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))
# plot!(subplot=4, mass2, label="$(title2)", xlabel="Time Step", title="Centered Mass", titlefont=font(10), legend=:topright, ylims=(-1e-5, +1e-5))

# # Comparing energy
# energy1 = readdlm("$(file1)_energy.csv", ',')
# energy2 = readdlm("$(file2)_energy.csv", ',')
# plot!(subplot=5, energy1, label="$(title1)", xlabel="Time Step", title="Energy", titlefont=font(10), legend=:topright)
# plot!(subplot=5, energy2, label="$(title2)", xlabel="Time Step", title="Energy", titlefont=font(10), legend=:topright)

# L2 norm over time

