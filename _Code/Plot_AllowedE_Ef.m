%% Plot the electron and hole energies along with the conduction and valence bands and the fermi level
figure;
hold on;
plot(x,Ec, '-b');
plot(x,Ev, '-r');

for i = 1 : numel(E_allowed_holes)
   
    plot([x(1), x(end)], [E_allowed_holes(i), E_allowed_holes(i)], '-r');
    
end

for i = 1 : numel(E_allowed_electrons)
   
    plot([x(1), x(end)], [E_allowed_electrons(i), E_allowed_electrons(i)], '-b');
    
end

plot([x(1), x(end)], [Ef_, Ef_], '-g', 'LineWidth', 1);

xlabel('m');
ylabel('eV');
