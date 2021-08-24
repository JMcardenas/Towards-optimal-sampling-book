%-------------------------------------------------------------------------%
% Filename: fig_4_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 4
%-------------------------------------------------------------------------%
% Inputs:
% row_num - row number (either 1 or 2)
% col_num - column number (either 1, 2 or 3)
% dim - dimension
%-------------------------------------------------------------------------%

function fig_4_data(row_num,col_num,dim)

%%% Define main parameters %%%
space     = ' ';
ir_domain      = 1;          % ir_domain: 1: Hypercube| 2: Anullus | 3: Linear const.
err_grid_ratio = 20;         % ratio between maximum index set size and error grid size
num_trials     = 1;          % number of random trials

% figure's name
fig_name = ['fig_4','_',num2str(row_num),'_',num2str(col_num),'_',num2str(dim)];

% index set
if col_num == 1
    % total product
    index_type = 'TP';
    
elseif col_num == 2
    % total degree
    index_type = 'TD';
    
else
    % hyperbolic cross
    index_type = 'HC';
end

% Set the dimension, maximum desired number of basis functions and the values of m
N_max = 6000;

if dim == 1
    d = 1;
elseif dim == 2
    d = 2;
elseif dim == 3
    d = 4;
elseif dim == 4
    d = 8;
elseif dim == 5
    d = 16;
else
    error('invalid subfigure selected');
end

%%% Construcft index set and define N %%%

n = find_order(index_type,d,N_max);      % find the maximum polynomial order for the given index_type
I = generate_index_set(index_type,d,n);  % compute index set
N = size(I,2);                           % number of basis functions

% n values
n_values = round(linspace(1,n,10));
num_n    = length(n_values);
N_values = zeros(num_n,1);

% arrays for storing the errors
num_cases        = 1;                   % legendre and/or orthogonal
Theta_values     = zeros(num_n,num_trials,num_cases);
Big_Theta_values = zeros(num_n,num_trials,num_cases);

%%% Construct error grid %%%

if ir_domain == 1                       % sampling points for the error grid
    err_samp_type = 'uniform';
else
    err_samp_type = sprintf('uniform_ir_%d',ir_domain);
end

M             = ceil(err_grid_ratio*N);
err_grid      = generate_sampling_grid(err_samp_type,d,M);

%%% Main loop %%%

for l = 1:num_cases
    
    % Legendre
    if l == 1
        basis_type     = 'legendre';        % functions used in the measurement matrix A
        err_basis_type = basis_type;        % polynomials used in computing the error
        
        % QR-Legendre
    else
        basis_type     = 'orthogonal';
        err_basis_type = basis_type;
        I              = reorder_index_set(I);          % Reorder index set
    end
    
    
    for i = 1:num_n
        
        n = n_values(i);
        
        % Construct index set
        I = generate_index_set(index_type,d,n);  % compute index set
        N = size(I,2);                           % number of basis functions
        
        N_values(i,1) = N;
        
        
        % Construct error matrix and error vector
        if row_num == 2
            
            if isequal(basis_type,'legendre')
                A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);
            else
                A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);
                [A_err_grid,R] = qr(A_err_grid,0);
            end
            
            % Compute theta
            Theta     = (const_theta(A_err_grid*sqrt(M)))^2;
            
            disp(['Figure 4_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),...
                space,'polynomials = ',basis_type,space,'Index type = ',index_type,space,'N = ',num2str(N),...
                space,'Theta = ',num2str(Theta)]);
            
            % save data
            Theta_values(i,:,l)     = Theta;
            
        else
            % Compute Theta
            Big_Theta = (max(generate_intrinsic_weights('legendre',I)))^2;
            
            disp(['Figure 4_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),...
                space,'polynomials = ',basis_type,space,'Index type = ',index_type,space,'N = ',num2str(N),...
                space,'Big Theta = ',num2str(Big_Theta)]);
            
            % save data
            Big_Theta_values(i,:,l) = Big_Theta;
        end
    end
end

Theta_data     = zeros(num_n,num_trials,num_cases);
Big_Theta_data = zeros(num_n,num_trials,num_cases);

for l = 1:num_cases
    if row_num == 2
        Theta_data(:,:,l)     = Theta_values(:,:,l);
    else
        Big_Theta_data(:,:,l) = Big_Theta_values(:,:,l);
    end
end

%%% save data %%%

clear A_err_grid R opt_w b_err_grid err_grid

save(['../../data/Figure 4/',fig_name])

