% PS
close all;
% clear all;
clc;
format long;

% Execute Potential_Guess - the potential function guess

%% Set up
iteration_counter = 1;
Initialize % constants and simulation parameters

%% Solve Shrodinger's equation for hole and electron wavefunctions
PS_SolveShrodinger_Electrons
PS_SolveShrodinger_Holes

%% Find Fermi Level
Find_Ef; % Ef_
% Fermi level conditions...
CalcChargeDensity;

%% Solve Poisson's Equation to find the outputted conduction band
Poisson_Equ

%%   
nss_s(:,1) = ns;
pss_s(:,1) = ps;
fermi_energies(1) = Ef_;
charge_density_mat(:,1) = Charge_Density;
electron_eigens{1} = E_allowed_electrons;
hole_eigens{1} = E_allowed_holes;
Ecs_In(:,1) = Ec;
Ec_Out = NewEc; % the Ec that was made
Ecs_Out(:,1) = Ec_Out;
Errors(1) = (sum( ( abs( Ecs_In(:, 1) - Ecs_Out(:, 1) ) ) ) );
Electrostatic_Potentials(:, 1) = Electrostatic_Potential;
Electric_Fields(:, 1) = Electric_Field;

save(FILE_NAME);

iteration_counter = iteration_counter + 1;

if(diag_disp)
    Ef_
    max(ns)
end

%% Start self-consistent solution finding
% Main Loop
while( (sum( ( abs( Ecs_In(:, iteration_counter - 1) - Ecs_Out(:, iteration_counter - 1) ) ) ) ) / numel(x) > converge_criterion)
    
    if(iteration_counter > MAX_SELF_CONS_ITER)
       break; 
    end
    
	f_factor = f_factor_vec(sum(iteration_counter > loop_counter_vec) + 1);
        
	fprintf('Current Error: %4.8f\n', (sum( ( abs( Ecs_In(:, iteration_counter - 1) - Ecs_Out(:, iteration_counter - 1) ) ) ) ) / numel(x) );
	pause(1);
        
	% Create new band profile using a mix of the input and output bands
	Ec_In = (1 - f_factor) * Ecs_In(:,iteration_counter - 1) + f_factor * Ecs_Out(:,iteration_counter - 1); % The input Ec
	Ecs_In(:, iteration_counter) = Ec_In;
        
	% Set the Ec and Ev
	Ec = Ec_In;
	Ev = Ec - Eg_vec;
                
	PS_SolveShrodinger_Electrons
	PS_SolveShrodinger_Holes
	Find_Ef
	CalcChargeDensity
	Poisson_Equ
        
    if(diag_disp)
        if(mod(iteration_counter, loop_graph_counter) == 0)
            close all;
            figure; plot(x,Ec); hold on; plot(x,Ev);
        end
    end
    
	Ec_Out = NewEc; % the Ec that was made
	
    Ecs_Out(:, iteration_counter) = Ec_Out;  
	nss_s(:,iteration_counter) = ns;
	pss_s(:,iteration_counter) = ps;
	fermi_energies(iteration_counter) = Ef_;
	charge_density_mat(:,iteration_counter) = Charge_Density;
	electron_eigens{iteration_counter} = E_allowed_electrons;
	hole_eigens{iteration_counter} = E_allowed_holes;   
    Electrostatic_Potentials(:, iteration_counter) = Electrostatic_Potential;
    Electric_Fields(:, iteration_counter) = Electric_Field;
    Errors(iteration_counter) = (sum( ( abs( Ecs_In(:, iteration_counter - 1) - Ecs_Out(:, iteration_counter - 1) ) ) ) );
	
    if(mod(iteration_counter, save_counter) == 0)
    	save(FILE_NAME);
    end  
            
	iteration_counter = iteration_counter + 1;
    
    if(diag_disp)
        Ef_
        max(ns)
    end
    
end    

save(FILE_NAME);