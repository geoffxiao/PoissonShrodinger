% Basis Set Expansion Solver 
% planck constant in eV * s
% Outputs: 
%   E_allowed (in eV)
%   wavefunctions_allowed

% Inputs:
mass_particle; % in kg / q units, q = 1.6e-19 C
N; % Number of basis functions to use, the more the merrier :)...but slower
potential_function; % in eV
z_axis; % in meters









% Start
%-------------------------------------------------------------------------%
if(numel(mass_particle) == 1)
    mass_particle = repmat(mass_particle, numel(z_axis), 1);
end

h = z_axis(2) - z_axis(1); % grid spacing
L = z_axis(end) - z_axis(1); % for basis function formation

if( iteration_counter <= 1 ) % No need to recalculate this stuff
    % Set up the basis functions, use the particle in a box basis set (sine
    % functions)
    Basis_Functions = zeros(numel(z_axis), N); 
    Laplaced = zeros(numel(z_axis), N);
    for i = 1 : N

        Basis_Functions(:,i) = sqrt(2/L) * sin(i * pi * (z_axis - z_axis(1)) / L);
        Laplaced(:,i) = -Basis_Functions(:,i) * (i * pi / L)^2 * (-h_bar^2 / 2 ) * (1 / mass_particle(1));
        
    end
end

if( sum(gradient(mass_particle)) ~= 0 )
    % Form the secular matrix
    % Hamiltonian =  (-h_bar^2 / 2) * ( d/dx ( 1 / mass(x) * d/dx) + V(x) )
    Hamiltonian_Operator = @(wavefunction, mass_particle, potential_function, h) (-h_bar^2 / 2 ) *...
        gradient(gradient(wavefunction, h) ./ mass_particle, h) + potential_function .* wavefunction;

    Secular_Matrix = zeros(N, N);
    for i = 1 : N % column

        for j = i : N % row

            % Hamiltonian acts on wavefunction
            % H(row, column) --> basis(row) * H * basis(column)
            Secular_Matrix(j, i) = trapz(z_axis, Basis_Functions(:, j) .*...
                Hamiltonian_Operator(Basis_Functions(:, i), mass_particle, potential_function, h));
            Secular_Matrix(i, j) = Secular_Matrix(j, i); % symmetric matrix

        end

    end

% Constant mass 
else
   
    Secular_Matrix = zeros(N, N);
    for i = 1 : N % column

        for j = i : N % row

            Secular_Matrix(j, i) = trapz(z_axis, Basis_Functions(:,j) .* (Laplaced(:,i) + potential_function .* Basis_Functions(:,i)));
            Secular_Matrix(i, j) = Secular_Matrix(j, i);
            
        end

    end
    
end

% coefficients to make the wavefunctions
% eigenenergies
[coeffs energies] = eig(Secular_Matrix);
wavefunctions_allowed = Basis_Functions * coeffs;
E_allowed = diag(energies);