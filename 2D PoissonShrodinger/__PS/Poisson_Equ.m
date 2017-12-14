% Solving the 2-D Poisson equation by the Finite Difference
...Method 
% Numerical scheme used is a second order central difference in space
...(5-point difference)

%%
%Specifying parameters
x_axis;
y_axis;
Charge_Density;
epsilon_rs;

%%
fprintf('Poisson Equation: ');

RHS = -Charge_Density / ( permittivity * epsilon_rs );

pn = zeros(nx,ny) + Inf; %Preallocating pn, next guess

%%
% Boundary conditions
p=zeros(nx, ny);
p(:,1)=0;
p(:,end)=0;
p(1,:)=0;                  
p(end,:)=0;

%%
i = 2 : nx - 1;
j = 2 : ny - 1;

%Explicit iterative scheme with C.D in space (5-point difference)
while( abs(max(max(pn - p))) > 1e-4 )

    pn = p; % current
    p(i,j)=( (dy^2*(pn(i+1,j)+pn(i-1,j)) )+( dx^2*(pn(i,j+1)+pn(i,j-1)) )-( RHS(i,j)*dx^2*dy*2) )/( 2*(dx^2+dy^2) );
    
    %Boundary conditions 
    p(:,1)=0;
    p(:,ny)=0;
    p(1,:)=0;                  
    p(nx,:)=0;

end

% Solution = p
calc_potential = -p; % negative b/cE = q * V, q < 0 for electrons
BandOffset; % NewEc, NewEv

fprintf('Complete...\n');