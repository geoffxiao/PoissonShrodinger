PHYSICAL_CONSTANTS

FILE_NAME = 'LAO_STO_14LAO_Passivated'; % Save to

format long;
diag_disp = 1; % 0/1 = don't/do display stuff while running

% MKS units, except energy in eV
Temperature = 300;
converge_criterion = 1e-6;
No_Holes = 0; % set to 1 to have no holes
No_Electrons = 0;

h = 1e-11; % discretization spacing
x = [0 : h : 153.33e-9]'; % Global variable for the x axis
interfaces = [381; 781];
for i = 3 : 27
   if(mod(i,2) == 0) 
       interfaces(i) = interfaces( i - 1 ) + 400;   
   else
       interfaces(i) = interfaces( i - 1 ) + 381;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material dependent...
Eg1 = 5.6; % band gap in eV
Eg2 = 3.2; % band gap in eV
Eg_params = [Eg1, Eg2];
Eg_params = repmat(Eg_params, 1, 14); 
% Left material is above right material by this much
dEc12 = 2.5; % 1 is above 2 by this much
dEc_params = [dEc12 -dEc12];
dEc_params = repmat(dEc_params, 1, 14);
dEc_params(end) = [];

epsilon_r1 = 27; % LAO
epsilon_r2 = 300; % STO
epsilon_r_params = [epsilon_r1, epsilon_r2];
epsilon_r_params = repmat(epsilon_r_params, 1, 14);

% Electron and hole effective masses
M_es = 0.4 * ones(28, 1) * M;
M_hs = 1.2 * ones(28, 1) * M;

% Set up the 2D density of states parameters, can simply set M_2D_es = M_es
M_2D_es = M_es;
M_2D_hs = M_hs;

LAO = LAO_Setup(10,0);
LAO_f = flipud(LAO);
STO = zeros(400, 1);
Other_Charges = repmat([LAO; STO], 14, 1);
Other_Charges = [Other_Charges; STO; STO; STO; STO; STO; STO; STO; STO; STO; STO; STO];
Other_Charges = Other_Charges + [0; ones(10,1) * ( 0.5 / ( 3.8e-10 * 3.905e-10 * 3.905e-10 * 0.25) ); zeros(numel(x) - numel(1:10) - 1, 1)]; % 1/2 unit cell LAO of H+

%%%% No need to change anything below this line... :)

%% M changes with position -> different sheet densities
M_2D_electrons = repmat(M_2D_es(1), interfaces(1), 1);
for i = 1 : numel(interfaces) - 1
   
    M_2D_electrons = [M_2D_electrons;...
        repmat(M_2D_es(i + 1), interfaces(i + 1) - interfaces(i), 1)];
    
end
M_2D_electrons = [M_2D_electrons;...
    repmat(M_2D_es(end), numel(x) - interfaces(end), 1)];


M_2D_holes = repmat(M_2D_hs(1), interfaces(1), 1);
for i = 1 : numel(interfaces) - 1
   
    M_2D_holes = [M_2D_holes;...
        repmat(M_2D_hs(i + 1), interfaces(i + 1) - interfaces(i), 1)];
    
end
M_2D_holes = [M_2D_holes;...
    repmat(M_2D_hs(end), numel(x) - interfaces(end), 1)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon_rs = repmat(epsilon_r1, interfaces(1), 1);
for i = 1 : numel(interfaces) - 1
   
    epsilon_rs = [epsilon_rs;...
        repmat(epsilon_r_params(i + 1), interfaces(i + 1) - interfaces(i), 1)];
    
end
epsilon_rs = [epsilon_rs;...
    repmat(epsilon_r_params(end), numel(x) - interfaces(end), 1)];

dEcs = -cumsum(dEc_params);
dEc_vec = zeros(interfaces(1), 1);
for i = 1 : numel(interfaces) - 1
   
    dEc_vec = [dEc_vec;...
        repmat(dEcs(i), interfaces(i + 1) - interfaces(i), 1)];
    
end
dEc_vec = [dEc_vec;...
    repmat(dEcs(end), numel(x) - interfaces(end), 1)];

Eg_vec = repmat(Eg_params(1), interfaces(1), 1);
for i = 1 : numel(interfaces) - 1
   
    Eg_vec = [Eg_vec;...
        repmat(Eg_params(i + 1), interfaces(i + 1) - interfaces(i), 1)];
    
end
Eg_vec = [Eg_vec;...
    repmat(Eg_params(end), numel(x) - interfaces(end), 1)];


%% Mass smoothing
% smoothing has little effect
smoothing_distance = 1;
smoothing_factor = 10000;% decrease to get more smoothing
M_e_vector = MassSmoothing(M_es, smoothing_distance, smoothing_factor, interfaces, x);
M_h_vector = MassSmoothing(M_hs, smoothing_distance, smoothing_factor, interfaces, x);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Potential Guess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess_init = [ zeros(interfaces(end), 1); ...
    ( [interfaces(end) + 1 : numel(x)]' - interfaces(end) ) * 10 / ( numel(x) - interfaces(end) ) ];

i = 0;
guess_init = guess_init + [([1 : interfaces(i + 1)]' - interfaces(i + 1)) * (-0.01);...
                           zeros(numel(x) - interfaces(i + 1),1)];
for i = 1 : 2 : 26
    guess_init = guess_init + [zeros(interfaces(i),1);...
                          ([interfaces(i) : interfaces(i + 1)]' - interfaces(i + 1)) * (1 + (60-i)/10) / 400;...
                           zeros(numel(x) - interfaces(i + 1) - 1,1)];
    
end

calc_potential = guess_init;
BandOffset % get Ec and Ev
Ec = NewEc;
Ev = NewEv;

% load from file
Ec = load('Ec.txt'); % text file guess
Ec = Ec + 0;
Ev = Ec - Eg_vec;
plot(x, Ec); hold on; plot(x, Ev);


%% Settings
% the charges are just number per m^3. NOT COLOUMBS/M^3!!!
DOPED = 0; % undoped, instead we're putting dopants as other charges...but can modify to include dopants later
Complete_Ionization = 1;
Donor_Activation_Pct = 1;

Donor_Dopants = []; % number per m^3
    Ed = NaN; % eV, donor state ionization energy relative to conduction band edge
Acceptor_Dopants = []; % number per m^3
    Ea = NaN; % eV, acceptor state ionization energy relative to valence band edge

% adaptive self consistency used for iterative solution
loop_counter_vec = [1000 2000 1000 Inf];
f_factor_vec = [0.005 0.03 0.01 0.01]; % the new input potential merging effect
MAX_SELF_CONS_ITER = 1000; % Max number of times to do the self-consistent solution finding
loop_graph_counter = 1; % plot Ec and Ev every ... times

% adaptive basis set expansion
% when iteration counter exceeds the value, advances to use the next N val for basis set expansion
% Larger N_adapt -> slower but better accuracy
N_adapt = [1000 2000 2500 4000 5000 Inf];
N_adapt_vals = [500 500 500 300 400 250];
save_counter = 10; % save every 10 times