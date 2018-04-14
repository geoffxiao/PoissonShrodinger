function out = LAO_Setup(unit_cells, sign)
    % sign: 1 = LaO on the x = 0 side, LaO = +
    % unit_cells = number of unit cells 
    
    % Unit cells of LAO construction
    
    h = 1e-11; % spacing
    % Unit cell size LAO = 3.8 A, Gu's paper
    c_axis = 3.8e-10;
    a_axis = 3.905e-10; % STO is cubic...?

    % LaO = +1, AlO2 = -1
    % Alternating layers, in a 0.5 * c * a^2 arrangement
    
    x_LAO = [0 : h : c_axis * unit_cells]';
    LAO_Charges = zeros(numel(x_LAO), 1);

    for i = 1 : unit_cells * 2

        for j = 1 : round(c_axis * 0.5 / h)

            LAO_Charges( round(c_axis * 0.5 * (i - 1) / h) + j ) = (-1)^(sign + 1) * 1 / (c_axis * a_axis^2 * 0.5);

        end

        sign = sign + 1;

    end
    LAO_Charges(1 : round(c_axis * 0.5 / h) : numel(x_LAO)) = 0;

    out = LAO_Charges;
    
end