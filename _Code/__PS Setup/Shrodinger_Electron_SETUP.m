E_limit = min([Ec(end), Ec(1)]); % top of condunction band well
potential_function = Ec;
z_axis = x;
mass_particle = M_e_vector; % divide by q to get energy in eV
N = N_adapt_vals(sum(iteration_counter > N_adapt) + 1);