%-------------------------------------------------------------------------%
% Filename: fig_3_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 3
%-------------------------------------------------------------------------%
% Inputs:
% row_num - row number (either 1 or 2)
% col_num - column number (either 1, 2 or 3)
%-------------------------------------------------------------------------%

function fig_3_data(row_num,col_num)

%%% Define main parameters %%%
space          = ' ';
index_type     = 'HC';                % use hyperbolic cross index set

err_grid_ratio = 30;                  % ratio between maximum index set size and error grid size
num_trials     = 50;                  % number of random trials
num_cases      = 1;                   % number of cases
num_m_type     = 5;                   % number of m type of values

fig_name = ['fig_3','_',num2str(row_num),'_',num2str(col_num)];

% Set the dimension, maximum desired number of basis functions and the values of m
d        = 2;
N_values = 3:100:1003;

%%% Construcft index set and define N %%%

n = find_order(index_type,d,max(N_values));     % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n);               % compute index set
N = size(I,2);                                  % number of basis functions

%%% set n values and N values

[n_values,N_values] = find_BigN_values(N_values,[],d,1);

num_m      = length(n_values);
N_vals_aux = [0 N_values];

% Reorder index set

I = reorder_index_set(I,n_values);

%%% Construct error grid %%%

if row_num == 1                           % sampling points for the error grid
    err_samp_type = 'uniform';
else
    err_samp_type = sprintf('uniform_ir_%d',row_num);
end

M        = ceil(err_grid_ratio*N);
rng('default'); rng(1);                 % set the seed
err_grid = generate_sampling_grid(err_samp_type,d,M);

%%% Construct error matrix and error vector %%%

err_basis_type = 'legendre';

A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);

[Q_err_grid,R_err_grid] = qr(A_err_grid,0);
A_err_grid = Q_err_grid;

% Compute optimal measure

mu = abs(A_err_grid.^2);

% arrays for storing the errors

sv_data = zeros(num_m,num_trials,num_cases,num_m_type);
C_data  = zeros(num_m,num_trials,num_cases,num_m_type);

%%% Main loop %%%
for m_type = 1:5
    
    if m_type == 1
        k_ratios = ones(1,num_m);
        m_values = k_ratios.*N_values;
        
    elseif m_type == 2
        k_ratios = 2*ones(1,num_m);
        m_values = k_ratios.*N_values;
        
    elseif m_type == 3
        k_ratios = 3*ones(1,num_m);
        m_values = k_ratios.*N_values;
        
    elseif m_type == 4
        k_ratios = round(log(N_values));
        m_values = k_ratios.*N_values;
        
    elseif m_type == 5
        k_ratios = 2*round(log(N_values));
        m_values = k_ratios.*N_values;
        
    else
        disp('wron m values');
    end
    
    if col_num < 3       %l = col_num  %1:num_cases-1
        
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
        
        parfor i = 1:num_m
            
            m = m_values(i);
            k = k_ratios(i);
            N = N_values(i);
            
            % temporary arrays to allow parfor
            sv_temp = zeros(1,num_trials);
            
            for t = 1:num_trials
                
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
                    disp('wrong weights');
                end
                
                % generate measurement matrix
                A = W*A_err_grid(I_set,1:N)/sqrt(m);
                
                sv_min = min(svd(A));
                C_min  = 1/sv_min;
                
                disp(['Figure 3','_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),space,...
                    'm type = ',num2str(m_type),space,'sampling = ',samp_type,space,'m = ',num2str(m),space,...
                    'trial = ',num2str(t),space,'iter = ',num2str(i),space,'sv min = ', num2str(sv_min),space,...
                    'C = ',num2str(C_min)]);
                
                sv_temp(1,t) = sv_min;
                
            end
            
            sv_data(i,:,1,m_type) = sv_temp;
            C_data(i,:,1,m_type)  = 1./sv_temp;
            
        end
    end
    
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
            sv_temp = zeros(1,num_m);
            
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
                
                sv_min = min(svd(A));
                C_min  = 1/sv_min;
                
                disp(['Figure 3','_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),space,...
                    'm type = ',num2str(m_type),space,'sampling = ',samp_type,space,'m = ',num2str(m),space,...
                    'trial = ',num2str(t),space,'iter = ',num2str(i),space,'sv min = ', num2str(sv_min),space,...
                    'C = ',num2str(C_min)]);
                
                sv_temp(1,i) = sv_min;
                
            end
            
            sv_data(:,t,1,m_type) = sv_temp;
            C_data(:,t,1,m_type)  = 1./sv_temp;
            
        end
        
    end
end

%%% Save data %%%
clear A_err_grid R_err_grid Q_err_grid b_err_grid y_grid I_set A b c W err_grid mu opt_mu

save(['../../data/Figure 3/',fig_name])


end