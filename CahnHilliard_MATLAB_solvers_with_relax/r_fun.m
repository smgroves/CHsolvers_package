function r = r_fun(phi,phi0,r0,b,hx,hy,C0,Beta,dt,gamma0)
    bphi0 = fft2(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);
    bphi  = fft2(b.*phi);
    bphi  = hx*hy*bphi(1,1);

    E1 = fft2(f(phi0,gamma0));
    r = r0 + 1/2*(bphi - bphi0) - Beta*dt*r0*(r0 - sqrt(E1(1,1)*hx*hy + C0));
end