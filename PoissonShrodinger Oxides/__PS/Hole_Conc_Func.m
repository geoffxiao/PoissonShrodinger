% Hole concentration function
% in m^-3

N_2D_holes = M_2D_holes .* k_b * Temperature / (3.14 * h_bar^2); 
Hole_Conc_Function = zeros(numel(x), 1);
ps = 0;

if(numel(E_allowed_holes) > 0)
    
    Hole_Conc_Function = zeros(numel(wavefunctions_allowed_holes(:,1)), 1);

    for i = 1 : numel(E_allowed_holes)
       
        ps_i = N_2D_holes .* int_FD(E_allowed_holes(i), Ef_, Temperature);
        Hole_Conc_Function = ps_i .* (wavefunctions_allowed_holes(:,i)).^2 + Hole_Conc_Function;
        ps = ps_i + ps;
        
    end

end

clear ps_i