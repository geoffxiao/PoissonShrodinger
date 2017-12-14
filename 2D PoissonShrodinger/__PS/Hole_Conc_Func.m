% Hole concentration function
% in m^-3

N_2D_holes = sqrt(2 * M_2D_holes .* k_b * Temperature / (pi* h_bar^2)); 
Hole_Conc_Function = zeros(nx, ny);
ps = 0;

if(numel(E_allowed_holes) > 0)
    
    for i = 1 : numel(E_allowed_holes)
       
        ps_i = N_2D_holes .* int_FD(E_allowed_holes(i), Ef_, Temperature);
        Hole_Conc_Function = ps_i .* (wavefunctions_allowed_holes(:,:,i)).^2 + Hole_Conc_Function;
        ps = ps_i + ps;
        
    end

end

clear ps_i