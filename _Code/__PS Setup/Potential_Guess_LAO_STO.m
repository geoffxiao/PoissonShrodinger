PHYSICAL_CONSTANTS

FILE_NAME = ''; % Save to

format long;
diag_disp = 1; % 0/1 = don't/do display stuff while running

% MKS units, except energy in eV
Temperature = 300;
converge_criterion = 5e-6;
No_Holes = 1; % set to 1 to have no holes
No_Electrons = 0;

h = 1e-11; % discretization spacing
x = [0 : h : 100e-9]'; % Global variable for the x axis
interfaces = [1000; 1960]; % index of the interfaces in the material

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material dependent...
Eg1 = 4.4; % band gap in eV
Eg2 = 3.1; % band gap in eV
Eg3 = 5.8; % band gap in eV
Eg_params = [Eg1, Eg2, Eg3];

% Left material is above right material by this much
dEc12 = 0.8; % 1 is above 2 by this much
dEc23 = -2.5; % 2 is above 3 by this much
dEc_params = [dEc12, dEc23];

epsilon_r1 = 6; % CSO
epsilon_r2 = 18; % BLSO
epsilon_r3 = 22; % LSAT
epsilon_r_params = [epsilon_r1, epsilon_r2, epsilon_r3];

% Other charges in the calculation
DonorAct = 0.65;
Other_Charges_params = [0, 7.18e26, 0, 0] * DonorAct;

% Electron and hole effective masses
M_es = [0.1; 0.34; 0.4] * M;
M_hs = [1.2; 1; 1.3] * M; % No holes so unimportant

% Set up the 2D density of states parameters, can simply set M_2D_es = M_es
M_2D_es = [M_es(1); M_es(2); M_es(3)];
M_2D_hs = [M_hs(1); M_hs(2); M_hs(3)];







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
calc_potential = zeros(numel(x), 1);
BandOffset % get Ec and Ev
Ec = NewEc;
Ev = NewEv;

% load from file
% Ec = load('Ec.txt'); % text file guess
% Ev = Ec - Eg_vec;
plot(x, Ec); hold on; plot(x, Ev);


%% Settings
% the charges are just number per m^3. NOT COLOUMBS/M^3!!!
DOPED = 0; % undoped, instead we're putting dopants as other charges...but can modify to include dopants later
Complete_Ionization = 1;
Donor_Activation_Pct = 1;

Other_Charges = repmat(Other_Charges_params(1), interfaces(1), 1);
for i = 1 : numel(interfaces) - 1
   
    Other_Charges = [Other_Charges;...
        repmat(Other_Charges_params(i + 1), interfaces(i + 1) - interfaces(i), 1)];
    
end
Other_Charges = [Other_Charges;...
    repmat(Other_Charges_params(end), numel(x) - interfaces(end), 1)];
Other_Charges = Other_Charges * Donor_Activation_Pct;

Donor_Dopants = []; % number per m^3
    Ed = NaN; % eV, donor state ionization energy relative to conduction band edge
Acceptor_Dopants = []; % number per m^3
    Ea = NaN; % eV, acceptor state ionization energy relative to valence band edge

% adaptive self consistency used for iterative solution
loop_counter_vec = [100 200 200 Inf];
f_factor_vec = [0.02 0.03 0.05 0.01]; % the new input potential merging effect
MAX_SELF_CONS_ITER = 200; % Max number of times to do the self-consistent solution finding
loop_graph_counter = 1; % plot Ec and Ev every ... times

% adaptive basis set expansion
% when iteration counter exceeds the value, advances to use the next N val for basis set expansion
% Larger N_adapt -> slower but better accuracy
N_adapt = [500 2000 2500 4000 5000 Inf];
N_adapt_vals = [15 500 200 300 400 250];
save_counter = 10; % save every 10 times