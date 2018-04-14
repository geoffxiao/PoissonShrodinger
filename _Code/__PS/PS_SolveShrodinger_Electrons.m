% Numerov Method implementation
% E_allowed_electrons
% wavefunctions_allowed_electrons

%% 
fprintf('\nShrodinger for Electrons: ');
    
Shrodinger_Electron_SETUP;

if(~No_Electrons)

    Shrodinger_BasisSetExpansion;
    
     % clean up, remove states with energy above the well
    if(~isempty(E_allowed))
        % clean up, remove states with energy above the well
        E_allowed = E_allowed(E_allowed < E_limit);
        wavefunctions_allowed = wavefunctions_allowed(:, 1 : numel(E_allowed));
        for i = 1 : numel(E_allowed)

            wavefunctions_allowed(:,i) = wavefunctions_allowed(:,i) / sqrt(trapz(z_axis, wavefunctions_allowed(:,i).^2));

        end 
    end

    wavefunctions_allowed_electrons = wavefunctions_allowed;
    E_allowed_electrons = E_allowed;
    
else
    
    E_allowed_electrons = [];
    wavefunctions_allowed_electrons = [];
    
end

fprintf('Complete \n');
