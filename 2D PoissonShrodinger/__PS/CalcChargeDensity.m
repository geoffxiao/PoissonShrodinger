% Calculate charge density

% units of p(z) are m^-3
Hole_Conc_Func % Hole_Conc_Function, ps

% units of n(z) are m^-3
Electron_Conc_Func % Electron_Conc_Function, ns

IonizedCharges; % calculate ionized charges

Charge_Density = q * (Hole_Conc_Function - Electron_Conc_Function + Other_Charges);  