% Poisson's Equation to find the potential
% Electric Field in V/m
% Electric Potential is in V
fprintf('Poisson Equation: ');

% del * [epsilon * del(potential)] = -charge density
% Electric Field = -del(potential)
% del(potential) = integral[ -charge density ] ./ epsilon
% Electric Field = integral[ charge density ] ./ epsilon
% potential = -integral(Electric Field)
% But we want to change potential to Ec, so we multiply by -1

Electric_Field = cumtrapz(z_axis, Charge_Density) ./ (permittivity * epsilon_rs);
Electrostatic_Potential = -cumtrapz(z_axis, Electric_Field);

% negative because E = q * V, and q < 0 for electrons
calc_potential = -Electrostatic_Potential;
BandOffset; % NewEc, NewEv

fprintf('Complete...\n');
