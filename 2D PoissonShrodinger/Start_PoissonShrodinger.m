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
Ecs_In(:,:,1) = Ec;
Ec_Out = NewEc; % the Ec that was made
Ecs_Out(:,:,1) = Ec_Out;
Errors(1) = (sum( ( abs( Ecs_In(:,:, 1) - Ecs_Out(:,:, 1) ) ) ) );

save(FILE_NAME);

iteration_counter = iteration_counter + 1;

if(diag_disp)
    Ef_
    max(ns)
end

%% Start self-consistent solution finding
% Main Loop
while( (sum( ( abs( Ecs_In(:,:, iteration_counter - 1) - Ecs_Out(:,:, iteration_counter - 1) ) ) ) ) / (nx * ny) > converge_criterion)
    
    if(iteration_counter > MAX_SELF_CONS_ITER)
       break; 
    end
    
	f_factor = f_factor_vec(sum(iteration_counter > loop_counter_vec) + 1);
        
	fprintf('Current Error: %4.8f\n', (sum( ( abs( Ecs_In(:,:, iteration_counter - 1) - Ecs_Out(:,:, iteration_counter - 1) ) ) ) ) / numel(x) );
	pause(1);
        
	% Create new band profile using a mix of the input and output bands
	Ec_In = (1 - f_factor) * Ecs_In(:,:,iteration_counter - 1) + f_factor * Ecs_Out(:,:,iteration_counter - 1); % The input Ec
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
        end
    end
    
	Ec_Out = NewEc; % the Ec that was made
	
    Ecs_Out(:,:, iteration_counter) = Ec_Out;  
    Errors(iteration_counter) = (sum( ( abs( Ecs_In(:,:, iteration_counter - 1) - Ecs_Out(:,:, iteration_counter - 1) ) ) ) );
	
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