%-------------------------------------------------------------------------%
% Filename: fig_1_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 1
%-------------------------------------------------------------------------%
% Inputs:
% row_num - row number (either 1, 2 or 3)
% col_num - column number (either 1, 2 or 3)
% dim     - dimension
%-------------------------------------------------------------------------%

function fig_1_data(row_num,col_num,dim)

%%% Define main parameters %%%
space          = ' ';
index_type     = 'HC';                % use hyperbolic cross index set

err_grid_ratio = 30;                  % ratio between maximum index set size and error grid size
num_trials     = 50;                   % number of random trials
num_cases      = 1;                   % number of cases
max_div        = 10;                  % number of m values

fig_name = ['fig_1','_',num2str(row_num),'_',num2str(col_num),'_',num2str(dim)];

% Set the function
if row_num == 1
    func      = @iso_exp;
    ir_domain = 2;
    
elseif row_num == 2
    func      = @split_product;
    ir_domain = 3;
    
elseif row_num == 3
    func      = @abs_rec;
    ir_domain = 2;
    
else
    error('invalid figure selected');
end

% Set the dimension, maximum desired number of basis functions and the values of m
N_max    = 1000;

if dim == 1
    d = 2;
    
elseif dim == 2
    d = 3;
    
elseif dim == 3
    d = 5;
    
elseif dim == 4
    d = 10;
    
elseif dim == 5
    d = 16;
    
else
    error('invalid subfigure selected');
end

%%% Construcft index set and define N %%%

n = find_order(index_type,d,N_max);     % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n);       % compute index set
N = size(I,2);                          % number of basis functions

%%% set m-values

n_values = find_n_values(N,d);
[avoid,N_values] = find_BigN_values([],n_values,d,2);
[m_values,avoid,N_values,k_ratios] = find_m_vals(d,n,max_div,N_values,2);

num_m      = length(m_values);
N_vals_aux = [0 N_values];

% Reorder index set

I = reorder_index_set(I,n_values);

%%% Construct error grid %%%

if ir_domain == 1                           % sampling points for the error grid
    err_samp_type = 'uniform';
else
    err_samp_type = sprintf('uniform_ir_%d',ir_domain);
end

M        = ceil(err_grid_ratio*N);
rng('default'); rng(1);                 % set the seed
err_grid = generate_sampling_grid(err_samp_type,d,M);

%%% Construct error matrix and error vector %%%

err_basis_type = 'legendre';

A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);

[Q_err_grid,R_err_grid] = qr(A_err_grid,0);
A_err_grid = Q_err_grid;

b_err_grid = func(err_grid)/sqrt(M);

% Compute optimal measure

mu = abs(A_err_grid.^2);

% arrays for storing the errors

L2_error_data   = zeros(num_m,num_trials,num_cases);
Linf_error_data = zeros(num_m,num_trials,num_cases);

%%% Main loop %%%

if  col_num < 3   %1:num_cases-1
    
    l = col_num;
    
    % choose sampling
    if l == 1
        % uniform
        samp_type = 'uniform_grid';
        
    elseif l == 2
        % method 1: non-adaptive
        samp_type = 'nonadaptive';
        
    else
        disp('wrong sampling');
    end
    
    for i = 1:num_m
        
        m = m_values(i);
        k = k_ratios(i);
        N = N_values(i);
        
        % temporary arrays to allow parfor
        L2_error_temp   = zeros(1,num_trials);
        Linf_error_temp = zeros(1,num_trials);
        
        parfor t = 1:num_trials
            
            % sample pts
            if isequal(samp_type,'uniform_grid')
                [y_grid, I_set] = generate_sampling_grid(samp_type,d,m,[],err_grid);
                
            elseif isequal(samp_type,'nonadaptive')
                opt_mu_non_adap = sum(mu(:,1:N),2)/N;
                [y_grid, I_set] = generate_sampling_grid(samp_type,d,m,opt_mu_non_adap,err_grid);
                
            else
                disp('wrong sampling');
            end
            
            % generate weights
            if isequal(samp_type,'uniform_grid')
                W = eye(m)*sqrt(M);
                
            elseif isequal(samp_type,'nonadaptive')
                W = 1./opt_mu_non_adap(I_set);
                W = diag(sqrt(W));
                
            else
                disp(['wrong weights']);
            end
            
            % generate measurement matrix
            A = W*A_err_grid(I_set,1:N)/sqrt(m);
            
            % generate measurement vector
            b = W*b_err_grid(I_set)/sqrt(m);
            
            % Solve system and save min sv
            c      = A\b;
            %sv_min = min(svd(A));
            
            % compute L^2_rho-norm
            
            L2_err  = norm(A_err_grid(:,1:N)*c - b_err_grid)/norm(b_err_grid);
            Linf_err =  norm(A_err_grid(:,1:N)*c - b_err_grid,Inf)/norm(b_err_grid,Inf);
            
            disp(['Figure 1','_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),space,...
                'sampling = ',samp_type,space,'m = ',num2str(m),space,'trial = ',num2str(t),space,'iter = ',...
                num2str(i),space,'L2 error = ', num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
            
            L2_error_temp(1,t)   = L2_err;
            Linf_error_temp(1,t) = Linf_err;
            
        end
        
        L2_error_data(i,:,1)   = L2_error_temp;
        Linf_error_data(i,:,1) = Linf_error_temp;
        
    end
end

%for l = num_cases:num_cases
if col_num == 3
    
    l = col_num;
    
    if l == 3
        % method 2: adaptive
        samp_type = 'adaptive';
        opt_mu    = mu;
        
    else
        disp('wrong sampling');
    end
    
    
    for t = 1:num_trials
        
        I_set = [];
        
        % temporary arrays to allow parfor
        L2_error_temp   = zeros(1,num_m);
        Linf_error_temp = zeros(1,num_m);
        
        for i = 1:num_m
            
            m = m_values(i);
            k = k_ratios(i);
            N = N_values(i);
            
            % sample pts
            if isequal(samp_type,'adaptive')
                for j = 1:N
                    
                    if j <= N_vals_aux(i)
                        k_sam = k_ratios(i)-k_ratios(i-1);
                    else
                        k_sam = k;
                    end
                    
                    [y_grid, I_new] = generate_sampling_grid(samp_type,d,k_sam,opt_mu(:,j),err_grid);
                    I_set           = [I_set; I_new];
                end
                
            else
                disp('wrong sampling');
            end
            
            % generate weights
            if isequal(samp_type,'adaptive')
                W = 1./(sum(mu(I_set,1:N),2)/N);
                W = diag(sqrt(W));
            else
                disp('wrong weights');
            end
            
            % generate measurement matrix
            A = W*A_err_grid(I_set,1:N)/sqrt(m);
            
            % generate measurement vector
            b = W*b_err_grid(I_set)/sqrt(m);
            
            % Solve system and save min sv
            c      = A\b;
            %sv_min = min(svd(A));
            
            % compute L^2_rho-norm
            
            L2_err   = norm(A_err_grid(:,1:N)*c - b_err_grid)/norm(b_err_grid);
            Linf_err =  norm(A_err_grid(:,1:N)*c - b_err_grid,Inf)/norm(b_err_grid,Inf);
            
            disp(['Figure 1','_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),space,...
                'sampling = ',samp_type,space,'m = ',num2str(m),space,'trial = ',num2str(t),space,'iter = ',...
                num2str(i),space,'L2 error = ', num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
            
            L2_error_temp(1,i)   = L2_err;
            Linf_error_temp(1,i) = Linf_err;
            
        end
        
        L2_error_data(:,t,1)   = L2_error_temp;
        Linf_error_data(:,t,1) = Linf_error_temp;
        
    end
    
end

%%% Save data %%%
clear A_err_grid Q_err_grid R_err_grid b_err_grid y_grid I_set A b c W err_grid mu opt_mu

save(['../../data/Figure 1/',fig_name])


end