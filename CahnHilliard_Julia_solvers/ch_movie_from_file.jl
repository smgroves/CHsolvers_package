# import Pkg;
# Pkg.add("ArgParse");
# Pkg.add("ColorSchemes");

using ArgParse
using Plots
using ColorSchemes
using VideoIO
using CSV
using DataFrames

function ch_movie_from_file(phi_file, t_out, ny; dtframes=1, filename="ch_movie", filetype="mp4", colorbar_type="default")
    """
    This function creates a red-white-blue video trajectory of chemical states.

    # Arguments
    - phi_file::String: File name for multidimensional array of chemical states.
    - t_out::Vector{Float64}: Vector of time steps corresponding to the third dimension of phi_t.
    - ny::Int: Number of mesh points in the y direction.

    # Optional Keyword Arguments
    - dtframes::Int: Number of dt time steps per frame. Default is 1.
    - filename::String: Name of the movie file to be saved. Default is 'ch_movie'.
    - filetype::String: Format of the movie file ('mp4', 'avi'). Default is 'mp4'.
    - colorbar_type::String: Type of colorbar to be used ('default', 'initial_range'). Default is 'default'.

    # Output
    Saves a red-white-blue video.
    """

    # Validate inputs
    valid_filetypes = ["mp4", "avi"]
    if !(filetype in valid_filetypes)
        error("Invalid filetype. Supported types are: mp4, avi")
    end

    # Open the video writer
    video_writer = VideoIO.open("$filename.$filetype")

    # Initial color range
    initial_range = nothing

    # Loop over time steps
    for i in 1:dtframes:length(t_out)
        # Read the necessary frame from the file
        row_start = (i - 1) * ny + 1
        row_end = i * ny
        phi_temp = CSV.read(phi_file, DataFrame; header=false, limit=row_end - row_start + 1, skip=row_start - 1) |> Matrix

        # Create the plot
        p = heatmap(
            phi_temp,
            color=:RdBu,
            colorbar=false,
            clims=begin
                if colorbar_type == "default"
                    (-1, 1)
                elseif colorbar_type == "initial_range"
                    initial_range === nothing && (initial_range = extrema(phi_temp))
                    initial_range
                elseif colorbar_type == "variable"
                    extrema(phi_temp)
                else
                    error("Invalid colorbar_type: $colorbar_type")
                end
            end,
            xlabel="",
            ylabel="",
            title="t = $(t_out[i])",
            aspect_ratio=1
        )

        # Write the frame to the video
        frame = Plots.frame(p)
        VideoIO.write(video_writer, frame[:color])

        # Print progress
        if mod(i, div(length(t_out), 20)) == 0
            println("$(round(i / length(t_out) * 100, digits=1))% complete")
        end
    end

    # Close the video writer
    close(video_writer)
end

# Example usage (adjust file paths and parameters as needed)
# ch_movie_from_file("phi_data.csv", 0:0.1:10, 128; dtframes=1, filename="output", filetype="mp4", colorbar_type="default")


# using Plots
# gr()
# x = 0:0.01:2*pi

# anim = @animate for i in 1:200
#     plot(x, sin.(x .+ i / 10.0), legend=false, size=(400, 200))
# end every 2;
# # gif(anim, "anim.gif") # then you can also change to  .mp4 or .mov
# #gif(anim,"anim.mp4") # this one usually works. If not, 
# #then you need to check ffmpeg. 

# using VideoIO
# using Plots

# @userplot CirclePlot
# @recipe function f(cp::CirclePlot)
#     x, y, i = cp.args
#     n = length(x)
#     inds = circshift(1:n, 1 - i)
#     linewidth --> range(0, 10, length=n)
#     seriesalpha --> range(0, 1, length=n)
#     aspect_ratio --> 1
#     label --> false
#     x[inds], y[inds]
# end

# n = 150
# t = range(0, 2π, length=n)
# x = sin.(t)
# y = cos.(t)

# anim = @animate for i ∈ 1:n
#     circleplot(x, y, i)
# end
# gif(anim, "anim_fps15.gif", fps=15)