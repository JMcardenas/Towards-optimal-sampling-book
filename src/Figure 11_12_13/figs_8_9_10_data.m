%-------------------------------------------------------------------------%
% Filename: figs_8_9_10_data.m 
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for sparse 
% approximation in high dimensions", Springer, 2021. 
%
% Description: generates and saves the data for Figures 8, 9 and 10.
%-------------------------------------------------------------------------%
% Inputs:
% fig_num - figure number (either 8, 9 or 10)
% row_num - row number (either 1 or 2)
% col_num - column number (either 1 or 2)
%-------------------------------------------------------------------------%

function figs_8_9_10_data(fig_num,row_num,col_num)

%%% Define main parameters %%%
space      = ' ';
index_type = 'HC';                    % use hyperbolic cross index set
 
err_grid_ratio = 10;                  % ratio between maximum index set size and error grid size
num_cases      = 2;                   % number of cases  

% figure's name
fig_name = ['fig','_',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

% Set the function
if fig_num == 8 
    func      = @iso_exp; 
    ir_domain = 2;
    
elseif fig_num == 9
    func      = @iso_exp; 
    ir_domain = 3;
    
elseif fig_num == 10 
    func      = @split_product; 
    ir_domain = 1;
    
else
    error('invalid figure selected');
end 

% Set the dimension, maximum desired number of basis functions and the values of m
 
if  col_num == 1
    d     = 2;
    N_max = 800;  
    
elseif col_num == 2
    d     = 8;
    N_max = 2000;  

elseif col_num == 3
    d     = 16;
    N_max = 4500;  
    
else
    error('invalid subfigure selected');
end
 

%%% Construcft index set and define N %%%

n = find_order(index_type,d,N_max);        % find the maximum polynomial order for the given index_type
I = generate_index_set(index_type,d,n);    % compute index set
N = size(I,2);                             % number of basis functions


%%% Construct error grid %%%

if ir_domain == 1                           % sampling points for the error grid
    err_samp_type = 'uniform';              
else
    err_samp_type = sprintf('uniform_ir_%d',ir_domain);      
end

M  = ceil(err_grid_ratio*N);

rng('default')
s = rng;

err_grid = generate_sampling_grid(err_samp_type,d,M);

% arrays for storing the errors
coeff_ref  = zeros(N,num_cases);  
coeff_sort = zeros(N,num_cases);
residual   = zeros(1,num_cases);


for l = 1:num_cases
    
    % Legendre, uniform
    if l == 1 
        basis_type     = 'legendre';        % functions used in the measurement matrix A
        err_basis_type = basis_type;        % polynomials used in computing the error
        
    % QR-Legendre, uniform
    else 
        basis_type     = 'orthogonal';
        err_basis_type = basis_type;
        
        % Reorder index set
        if row_num == 2
            I = reorder_index_set(I,[]);   
        end
    end
    
    %%% Construct error matrix and error vector %%%
    
    A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid); 
    
    if isequal(basis_type,'orthogonal') 
        [Q_err_grid,R] = qr(A_err_grid,0); 
        A_err_grid     = Q_err_grid;
    end
    
    b_err_grid = func(err_grid)/sqrt(M);    
    
    %%% Compute a reference solution via least squares on the error grid %%%
    
    coeff_ref(:,l) = A_err_grid\b_err_grid;
    residual(l)    = norm(A_err_grid*coeff_ref(:,l) - b_err_grid);
    
    %%% Sort coeff %%%
    
    [coeff_sort(:,l),J] = sort(abs(coeff_ref(:,l)),'descend');
    
    disp(['Figure ',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),...
             space,'dimension = ',num2str(d),space,'polynomials = ',basis_type,...
             space,'re order = ',num2str(row_num)]);
    
end


%%% Save data %%%
clear A_err_grid R b_err_grid Q_err_grid

save(['../../data/Figure 8_9_10/',fig_name])

end
