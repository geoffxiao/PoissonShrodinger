nss_s = zeros(numel(x), MAX_SELF_CONS_ITER);
pss_s = zeros(numel(x), MAX_SELF_CONS_ITER);

charge_density_mat = zeros(numel(x), MAX_SELF_CONS_ITER);

electron_eigens = cell(MAX_SELF_CONS_ITER,1);
hole_eigens = cell(MAX_SELF_CONS_ITER,1);

fermi_energies = zeros(MAX_SELF_CONS_ITER, 1);

Ecs_Out = zeros(numel(x), MAX_SELF_CONS_ITER);
Ecs_In = zeros(numel(x), MAX_SELF_CONS_ITER);

Errors = zeros(MAX_SELF_CONS_ITER, 1);

Electrostatic_Potentials = zeros(numel(x), MAX_SELF_CONS_ITER);
Electric_Fields = zeros(numel(x), MAX_SELF_CONS_ITER);

% z_axis and potential_function only pertain to the Shrodinger solver
% part...globally:
    % x = the x axis
    % Ec, Ev = the conduction and valence bands