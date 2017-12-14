[x, y] = meshgrid(x_axis, y_axis);

dx = x_axis(2) - x_axis(1); % grid spacing
L_x = x_axis(end) - x_axis(1); % for basis function formation

dy = y_axis(2) - y_axis(1); % grid spacing
L_y = y_axis(end) - y_axis(1); % for basis function formation

nx = numel(x_axis);
ny = numel(y_axis);

% Need to order the basis functions by energy!!
% N must be perfect square
Ordering = zeros(N,2); 
Ordering(:,1) = repmat([1:round(sqrt(N))],1,round(sqrt(N)))';
c = 1;
for i = 1 : round(sqrt(N))
    for j = 1 : round(sqrt(N))
        v(c) = i;
        c = c + 1;
    end
end
v = v';
Ordering(:,2) = v;
Ordering(:,3) = (Ordering(:,1).^2 / L_x^2) + (Ordering(:,2).^2 / L_y^2);
[values, order] = sort(Ordering(:,3));
QuantumNumbers = Ordering(order,:);

% Set up the basis functions, use the particle in a box basis set (sine
% functions)
Basis_Functions = zeros(numel(x_axis), numel(y_axis), N); 
Laplaced = zeros(numel(x_axis), numel(y_axis), N); 

% Construct basis functions and Laplacian operator on basis functions
N_loop = round(sqrt(N));
counter = 1;
for counter = 1 : N
        
    i = QuantumNumbers(counter,1);
    j = QuantumNumbers(counter,2);
    Basis_Functions(:,:,counter) = [( sqrt(2/L_x) * sin((i * pi * (x - x_axis(1))) / L_x) ) .* ...
                                           ( sqrt(2/L_y) * sin((j * pi * (y - y_axis(1))) / L_y) )]';
    Laplaced(:,:,counter) = -Basis_Functions(:,:,counter) * ( (i * pi / L_x)^2 + (j * pi / L_y)^2 ) * (-(h_bar^2) / (2 * mass_particle) ); 
    
end

%%
Ecs_Out = zeros(numel(x_axis), numel(y_axis), MAX_SELF_CONS_ITER);
Ecs_In = zeros(numel(x_axis), numel(y_axis), MAX_SELF_CONS_ITER);

Errors = zeros(MAX_SELF_CONS_ITER, 1);