function err = compute_L2(phi, phi_ref)
% Grid-based L2 error (no h² factor, matching group convention).
diff = phi - phi_ref;
err  = sqrt(sum(diff(:).^2));
end