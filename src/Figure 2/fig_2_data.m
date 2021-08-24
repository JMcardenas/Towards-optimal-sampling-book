%-------------------------------------------------------------------------%
% Filename: fig_2_data.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 2
%-------------------------------------------------------------------------%
% Inputs:
% row_num - row number (either 1 or 2)
% col_num - column number (either 1, 2 or 3)
%-------------------------------------------------------------------------%

function fig_2_data(row_num,col_num)

%%% Define main parameters %%%
space          = ' ';
index_type     = 'HC';                % use hyperbolic cross index set

err_grid_ratio = 20;                  % ratio between maximum index set size and error grid size
num_trials     = 1;                   % number of random trials
num_cases      = 2;                   % number of cases
max_div        = 20;                  % number of m values

fig_name = ['fig_2','_',num2str(row_num),'_',num2str(col_num)];


% Set the dimension, maximum desired number of basis functions and the values of m
N_max = 1000;
d     = 2;

%%% Construcft index set and define N %%%

n = find_order(index_type,d,N_max);     % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n);       % compute index set
N = size(I,2);                          % number of basis functions

%%% set m-values

[m_values,n_values,N_values,k_ratios] = find_m_vals(d,n,max_div,[],1);
num_m  = length(m_values);

% Reorder index set

I = reorder_index_set(I,n_values);

%%% Construct error grid %%%

if col_num == 1                           % sampling points for the error grid
    err_samp_type = 'uniform';
else
    err_samp_type = sprintf('uniform_ir_%d',col_num);
end

M        = ceil(err_grid_ratio*N_max);
err_grid = generate_sampling_grid(err_samp_type,d,M);

%%% Construct error matrix and error vector %%%

err_basis_type = 'legendre';

A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);

[Q_err_grid,R_err_grid] = qr(A_err_grid,0);
A_err_grid = Q_err_grid;

% Compute optimal measure

mu = abs(A_err_grid.^2);

% arrays for storing the errors

I_set_data = [];

%%% Main loop %%%

for l = row_num
    
    % choose sampling
    if l == 1
        % uniform
        samp_type = 'uniform_grid';
        
    elseif l == 2
        % method 1: non-adaptive
        samp_type = 'nonadaptive';
        
    elseif l == 3
        % method 2: adaptive
        samp_type = 'adaptive';
        opt_mu    = mu;
        
    else
        disp('wrong sampling');
    end
    
    I_set = [];
    
    for i = 1:8
        
        m = m_values(i);
        k = k_ratios(i);
        N = N_values(i);
        
        
        % sample pts
        if isequal(samp_type,'uniform_grid')
            [y_grid, I_set] = generate_sampling_grid(samp_type,d,m,[],err_grid);
            
        elseif isequal(samp_type,'nonadaptive')
            opt_mu_non_adap = sum(mu(:,1:N)')/N;
            [y_grid, I_set] = generate_sampling_grid(samp_type,d,m,opt_mu_non_adap,err_grid);
            
        elseif isequal(samp_type,'adaptive')
            for j = 1:N
                
                if (i > 1)
                    if j <= N_values(i-1)
                        k_sam = k_ratios(i)-k_ratios(i-1);
                    end
                else
                    k_sam = k;
                end
                
                [y_grid, I_new] = generate_sampling_grid(samp_type,d,k_sam,opt_mu(:,j),err_grid);
                I_set           = [I_set; I_new];
            end
            
        else
            disp('wrong sampling');
        end
        
        disp(['Figure 2','_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),...
            space,'sampling = ',samp_type,space,'m = ',num2str(m),space,'iter = ',num2str(i)]);
        
        %I_set_data = [I_set_data; I_set];
        
    end
end


%%% Save data %%%
clear A_err_grid R_err_grid y_grid Q_err_grid mu


save(['../../data/Figure 2/',fig_name])


end
