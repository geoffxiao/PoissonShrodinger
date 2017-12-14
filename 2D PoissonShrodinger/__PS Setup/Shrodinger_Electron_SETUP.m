E_limit = min([Ec(end,:), Ec(1,:), Ec(:,1), Ec(:,end)]); % top of condunction band well
potential_function = Ec;

mass_particle = M_e_vector; % divide by q to get energy in eV