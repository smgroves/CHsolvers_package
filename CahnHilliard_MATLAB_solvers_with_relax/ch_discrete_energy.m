function E = ch_discrete_energy(phi,hxy,eps2,gamma0)
    [gridx,gridy] = size(phi);
    a = hxy*sum(sum(f(phi,0))); %Calculate chemical free energy
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
    end
    E = a + 0.5*eps2*(sum_i+sum_j);
end
% function E = ch_discrete_energy(phi, hxy, eps2, gamma0)
%     % Get grid dimensions
%     [gridx, gridy] = size(phi);
    
%     % Calculate chemical free energy (unchanged)
%     a = hxy * sum(sum(f(phi, 0)));
    
%     % Compute discrete Laplacian with periodic boundaries
%     phi_right = circshift(phi, [0, 1]);   % Shift right (y-direction)
%     phi_left = circshift(phi, [0, -1]);   % Shift left (y-direction)
%     phi_up = circshift(phi, [-1, 0]);     % Shift up (x-direction)
%     phi_down = circshift(phi, [1, 0]);    % Shift down (x-direction)
%     Lap_phi = phi_right + phi_left + phi_up + phi_down - 4 * phi;
%     % Note: Assuming h=1; if h≠1, divide Lap_phi by h^2 and adjust hxy accordingly
    
%     % Compute interfacial free energy using integration by parts
%     interfacial_energy = - (eps2 / 2) * hxy * sum(sum(phi .* Lap_phi));
    
%     % Total energy
%     E = a + interfacial_energy;
% end