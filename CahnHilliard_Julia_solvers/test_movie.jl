# using Plots
# gr()
# x = 0:0.01:2*pi

# anim = @animate for i in 1:200
#     plot(x, sin.(x .+ i / 10.0), legend=false, size=(400, 200))
# end every 2;
# # gif(anim, "anim.gif") # then you can also change to  .mp4 or .mov
# #gif(anim,"anim.mp4") # this one usually works. If not, 
# #then you need to check ffmpeg. 

using VideoIO
using Plots

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length=n)
    seriesalpha --> range(0, 1, length=n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

n = 150
t = range(0, 2π, length=n)
x = sin.(t)
y = cos.(t)

anim = @animate for i ∈ 1:n
    circleplot(x, y, i)
end
gif(anim, "anim_fps15.gif", fps=15)