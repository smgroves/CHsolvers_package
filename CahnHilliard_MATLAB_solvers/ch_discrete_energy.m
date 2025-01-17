function E = ch_discrete_energy(phi,hxy,gridx,gridy,eps2)
%Local function for calculating the discrete energy across the domain
% f = @(x) 0.25*(x.^2-1).^2; %Define double-well free energy
a = hxy*sum(sum(f(phi))); %Calculate chemical free energy
sum_i = 0; %Initialize interfacial free energy in x
for i = 1:gridx-1
    for j = 1:gridy
        sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
    end
end
sum_j = 0; %Initialize interfacial free energy in y
for i = 1:gridx
    for j = 1:gridy-1
        sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
    end
end
E = a + 0.5*eps2*(sum_i+sum_j);
end