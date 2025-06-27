include("../CahnHilliard_Julia_solvers/CahnHilliard_NMG.jl")
function initialize_round_CPC_droplet_only_um(nx, ny; CPC_width=0.173, domain_width=3.2, epsilon=0.1)
    CPC_radius_R0 = CPC_width / domain_width

    phi = zeros(Float64, nx, ny)
    h = 1 / nx
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    R = @.sqrt((xx - 0.5)^2 + (yy - 0.5)^2)
    # eps_c = epsilon
    # delta = eps_c * sqrt(2)
    # psi0 = 0.5 * (1 .+ @.tanh((R0 .- R) / (2 * delta)))
    # phi = 2 .* psi0 .- 1    # psi0=(phi0+1)/2
    phi = @.tanh((CPC_radius_R0 .- R) / (sqrt(2) * epsilon))
    return phi
end



epsilon = 0.0067
nx = 512
phi = initialize_round_CPC_droplet_only_um(nx, nx, CPC_width=0.8, domain_width=6.4, epsilon=0.0067)

writedlm("./IC/CPC_droplet_e_0.0067_128_R0_0.8.csv", phi, ',')
