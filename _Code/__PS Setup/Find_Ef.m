% Find Fermi Energy in eV
% Ef_ = Fermi Energy found
% ns, ps = electron, hole sheet density in m^-2
% m^-2 = (1/10000) cm^-2
% Overall charge neutrality
epsilon_Ef = 1e-30;
Ef_prev = Inf;
Ef_low = min(Ev) - max(Eg_vec) * 2;
Ef_high = max(Ec) + max(Eg_vec) * 2;

% Energy will be in eV
Ef_ = ( Ef_low + Ef_high ) / 2;

CalcChargeDensity;
Total_Charge = trapz(x, Charge_Density);
 
while (abs(Total_Charge) > epsilon_Ef) && (abs(Ef_prev - Ef_) > 1e-16)

	% too many electrons, Ef too high
	if(sign(Total_Charge) < 0)
        Ef_prev = Ef_;
        Ef_high = Ef_;
        Ef_ = (Ef_high + Ef_low) / 2;
    else
        Ef_prev = Ef_;
        Ef_low = Ef_;
        Ef_ = (Ef_high + Ef_low) / 2;
    end

	if(0)
        Ef_
    end 
        
	CalcChargeDensity;
	Total_Charge = trapz(x, Charge_Density);        
        
end

% CalcChargeDensity;
% Ef_;

clear Ef_guess