% Electron concentration function
% in m^-3
% ns

N_2D_electrons = M_2D_electrons .* k_b * Temperature / (3.14 * h_bar^2); % E_allowed_electrons
Electron_Conc_Function = zeros(numel(x), 1);
ns = 0;

if(numel(E_allowed_electrons) > 0)
    
    Electron_Conc_Function = zeros(numel(wavefunctions_allowed_electrons(:,1)), 1);

    for i = 1 : numel(E_allowed_electrons)

        ns_i = N_2D_electrons .* int_FD(Ef_, E_allowed_electrons(i), Temperature);
        Electron_Conc_Function = ns_i .* (wavefunctions_allowed_electrons(:,i)).^2 + Electron_Conc_Function;
        ns = ns_i + ns;
        
    end

end


clear ns_i