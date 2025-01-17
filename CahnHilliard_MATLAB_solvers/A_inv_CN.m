function Ai = A_inv_CN(phi,dt,k2,k4,gamma0,epsilon2)
    denom = 1 + dt/2*epsilon2*k4 - dt/2*gamma0*k2;
    Ai = real(ifft2(fft2_filtered(phi)./denom));
end