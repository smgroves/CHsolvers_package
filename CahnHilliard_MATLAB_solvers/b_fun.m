function b = b_fun(phi,hx,hy,C0,gamma0)
    E1 = fft2_filtered(f(phi));
    b = df(phi,gamma0) ./ sqrt(E1(1,1)*hx*hy + C0);
end