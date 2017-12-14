E_limit = min(-Ev(end), -Ev(1)); % bottom of the valence band well, negated
E_limit = min(-Ev(end)); % bottom of the valence band well, negated

potential_function = -Ev;
z_axis = x;
mass_particle = M_h_vector; % divide by q to get energy in eV
N = N_adapt_vals(sum(iteration_counter > N_adapt) + 1);