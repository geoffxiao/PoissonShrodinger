function [ out ] = int_FD(E, Ef, Temperature)

    k_b = 8.617e-5;
    
    % Fermi-Dirac
    
    
    % Prevent Inf overflow!
    if( isinf(  exp((E - Ef) / (k_b * Temperature)) ) )
        
        out = (E - Ef) / (k_b * Temperature) ; % log(1 + exp(x)) ~ log(exp(x)) ~ x for x-->inf
        
    % also if it is too small... 
    elseif(exp((E - Ef) / (k_b * Temperature)) < 1e-16)
        
        out = exp((E - Ef) / (k_b * Temperature)); % log(1 + x) = x for x-->0, x = exp(stuff) in this case
            
    else
            
        out = (log(   1 + exp(   (E - Ef) / (k_b*Temperature)  ) ));
        
    end

end
