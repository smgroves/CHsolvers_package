function r0 = r0_fun(phi0,hx,hy,C0)
    fphi = f(phi0);
    r0 = sqrt(hx*hy*sum(sum(fphi)) + C0);
end