% Basis Set Expansion Solver 
% planck constant in eV * s
% Outputs: 
%   E_allowed (in eV)
%   wavefunctions_allowed

% Inputs:
mass_particle; % in kg / q units, q = 1.6e-19 C
N = 15^2; % Number of basis functions to use, the more the merrier :)...but slower
x_axis; % in meters
y_axis; % in meters
potential_function;

%% Changes each time b/c potential function changes
% Hamiltonian operator operated on wavefunction
Hamiltonian = zeros(numel(x_axis), numel(y_axis), N); 

for i = 1 : N
    
    Hamiltonian(:,:,i) = Laplaced(:,:,i) + potential_function .* Basis_Functions(:,:,i);
    
end

%%
% Form the secular matrix
% Hamiltonian = (-h_bar^2 / 2*m) * Laplacian + V(x,y)
Secular_Matrix = zeros(N, N);
for i = 1 : N % column
    
    for j = i : N % row
   
        % Hamiltonian acts on wavefunction
        % H(row, column) --> basis(row) * H * basis(column)
        
        % Approximate 2D integral as sum of sum b/c each element represents one
        % grid point
        F = Basis_Functions(:,:,j) .* Hamiltonian(:,:,i);
        Secular_Matrix(j, i) = trapz(y_axis, trapz(x_axis, F));        
        Secular_Matrix(i, j) = Secular_Matrix(j, i); % symmetric matrix
        
    end
    
end

%%
% coefficients to make the wavefunctions
% eigenenergies
[coeffs, energies] = eig(Secular_Matrix); % eigs = 6, eig = all
wavefunctions_allowed = reshape(reshape(Basis_Functions, [], size(Basis_Functions, 3)) * coeffs, size(Basis_Functions));

E_allowed = diag(energies);

% wavefunctions_allowed = zeros(numel(x_axis), numel(y_axis), N_lim);
% for i = 1 : N_lim % calc all wavefunctions
%     for j = 1 : N
%         wavefunctions_allowed(:,:,i) = wavefunctions_allowed(:,:,i) + Basis_Functions(:,:,j) * coeffs(j,i);
%     end
% end