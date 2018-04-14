% Ionization of donors and acceptors

% Donors 
if(DOPED)
    if(~Complete_Ionization)
        Nd_Ionized = Donor_Dopants ./ (1 + exp( Ef_ - (Ec + Ed) ));
        Na_Ionized = Acceptor_Dopants ./ (1 + exp( (Ev + Ea) - Ef_) );
    else
        Nd_Ionized = Donor_Dopants;
        Na_Ionized = Acceptor_Dopants;
    end
    % acceptor states make material p-type, when they accept an electron the
    % state becomes negative b/c extra electron
    Ionized = Nd_Ionized - Na_Ionized;    
else
    Ionized = zeros(numel(x), 1);
end

