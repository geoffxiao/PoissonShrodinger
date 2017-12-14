E_limit = min(-Ev(end,:), -Ev(1,:), -Ev(:,1), -Ev(:,end)); % bottom of the valence band well, negated

potential_function = -Ev;
mass_particle = M_h_vector; % divide by q to get energy in eV