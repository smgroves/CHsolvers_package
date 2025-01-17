function r = r_fun(phi,phi0,r0,b,hx,hy,C0,Beta,dt)
    bphi0 = fft2_filtered(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);
    bphi  = fft2_filtered(b.*phi);
    bphi  = hx*hy*bphi(1,1);

    E1 = fft2_filtered(f(phi0));
    r = r0 + 1/2*(bphi - bphi0) - Beta*dt*r0*(r0 - sqrt(E1(1,1)*hx*hy + C0));
end