#%%
include("aux_functions_SAV.jl")

#%%
#make a 2 by 2 matrix x
x = [1.0 2.0; 3.0 4.0]

println(ext(x))

println(extback(ext(x)))

#%%
boundary = "neumann"
nx, ny = size(x)
xright, xleft, yright, yleft = [1 0 1 0]
Lx = xright - xleft
Ly = yright - yleft

# % Decide on the solver's mesh spacing for NEUMANN vs PERIODIC
# %  - For Neumann: we will mirror the domain, so pass 2*hx and 2*hy into sav_solver.
# %  - For Periodic: keep as-is.
if boundary == "neumann"
    Lx = 2 * Lx
    Ly = 2 * Ly
    nx = 2 * nx
    ny = 2 * ny
end
hx = Lx / nx
hy = Ly / ny
h2 = hx * hy #Define mesh size


#%%
m = 8
epsilon2 = NaN

if isnan(epsilon2)
    epsilon2 = h2 * m^2 / (2 * sqrt(2) * atanh(0.9))^2 #Define ϵ^2 if not prespecified
else
    m = sqrt((epsilon2 * (2 * sqrt(2) * atanh(0.9))^2) / h2) #Else overwrite m
end
#%%
kx = 1im * vcat(0:nx÷2, -nx÷2+1:-1) * (2π / Lx)
ky = 1im * vcat(0:ny÷2, -ny÷2+1:-1) * (2π / Ly)

kxx = kx .^ 2
kyy = ky .^ 2
meshgrid(x, y) = (repeat(x, 1, length(y)), repeat(y', length(x), 1))
kxx_mat, kyy_mat = meshgrid(kxx, kyy)

k2 = (kxx_mat + kyy_mat)
k4 = k2 .^ 2


#%%
r0_fun(x, hx, hy, 0)

#%%
df(x, 0)

f(x)

#%%
b_fun(x, hx, hy, 0, 0)
fft_filtered(x)


#%%
Lap_SAV(ext(x), k2)

#%%
dt = 1e-5
A_inv_CN(ext(x), dt, k2, k4, 0, epsilon2)

function A_inv_CN(phi, dt, k2, k4, gamma0, epsilon2)
    denom = 1 .+ (dt / 2) * epsilon2 .* k4 .- (dt / 2) * gamma0 .* k2
    return real.(ifft(fft_filtered(phi) ./ denom))
end


#%%
fft_filtered(x)

#%%
C0 = 0
r0 = r0_fun(ext(x), hx, hy, C0) # Initialize sav state
phi0_df = df(ext(x), gamma0) #df at phi0
Lap_dfphi0 = Lap_SAV(phi0_df, k2)    #Lap of df(phi0)
phi_bar = A_inv_CN(ext(x) + dt / 2 * Lap_dfphi0, dt, k2, k4, gamma0, epsilon2)

b = b_fun(phi_bar, hx, hy, C0, gamma0)
Beta = 0
g_fun_CN(ext(x), r0, b, dt, hx, hy, epsilon2, gamma0, Beta, C0, k2)
