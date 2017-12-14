% Numerov Method implementation
% For holes
% For the input potential, just put it upside down and solve the same way
% This works because band profile is energy 
% Energy = charge * potential
% charge of electron and hole are opposite!
% E_allowed_holes = []
% wavefunctions_allowed_holes
%%
% Ev
fprintf('Shrodinger for Holes: ');

Shrodinger_Hole_SETUP

if(~No_Holes)
    
    Shrodinger_BasisSetExpansion;
    
    % clean up, remove states with energy above the well
    if(~isempty(E_allowed))
        % clean up, remove states with energy above the well
        E_allowed = E_allowed(E_allowed < E_limit);
        wavefunctions_allowed = wavefunctions_allowed(:, 1 : numel(E_allowed)); % Already normalized
    end
    
    wavefunctions_allowed_holes = wavefunctions_allowed;
    E_allowed_holes = -E_allowed;
    
else
    
    E_allowed_holes = [];
    wavefunctions_allowed_holes = [];
    
end

fprintf('Complete \n');