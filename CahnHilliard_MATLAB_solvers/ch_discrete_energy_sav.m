function E = ch_discrete_energy_sav(phi,hxy,eps2,k2,gamma0,r,C0)
    % % % E1 = hxy*fft2(f(phi,0)); %Calculate chemical free energy
    % % E1 = hxy*fft2(f(phi,gamma0)); %Calculate chemical free energy
    % % E1 = E1(1,1);

    % E = 0.5 * fft2(phi .*(-eps2*real(ifft2(k2 .* fft2(phi)))));
    % E = r^2 + E(1,1);

    E_gamma = hxy*fft2(gamma0/2*phi.^2);
    E_gamma = E_gamma(1,1);

    [gridx,gridy] = size(phi);
    a = r^2-C0-(gamma0^2+2*gamma0)/4; 
    sum_i = 0; % Initialize interfacial free energy in x
    for i = 1:gridx-1
        for j = 1:gridy
            sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
        end
    end
    sum_j = 0; % Initialize interfacial free energy in y
    for i = 1:gridx
        for j = 1:gridy-1
            sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
    end
    E = a + 0.5*eps2*(sum_i+sum_j)+E_gamma;

end