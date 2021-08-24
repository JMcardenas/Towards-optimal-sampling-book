%-------------------------------------------------------------------------%
% Filename: figs_11_12_13_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 11, 12 and 13.
%-------------------------------------------------------------------------%
% Inputs:
% fig_num - figure number (either 11, 12 or 13)
% row_num - row number (either 1 or 2)
% col_num - column number (either 1 or 2)
%-------------------------------------------------------------------------%

function figs_11_12_13_data(fig_num,row_num,col_num)

%%% Define main parameters %%%
space          = ' ';
index_type     = 'HC';                % use hyperbolic cross index set

err_grid_ratio = 10;                  % ratio between maximum index set size and error grid size
num_trials     = 50;                  % number of random trials
num_cases      = 4;                   % number of cases
weighted_opt   = 0;                   % weighted_opt - 1: weighted | 0: non-weighted

% CVX options
cvx_opt.precision = 'default';
cvx_opt.solver    = 'mosek';
cvx_opt.verbose   = false;

fig_name = ['fig','_',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

% Set the function
if fig_num == 11
    func      = @iso_exp;
    ir_domain = 2;
    
elseif fig_num == 12
    func      = @iso_exp;
    ir_domain = 3;
    
elseif fig_num == 13
    func      = @split_product;
    ir_domain = 1;
    
else
    error('invalid figure selected');
end

% Set the dimension, maximum desired number of basis functions and the values of m
if row_num == 1 && col_num == 1
    d        = 1;
    N_max    = 400;
    m_values = 10:10:100;
    
elseif row_num == 1 && col_num == 2
    d        = 2;
    N_max    = 800;
    m_values = 20:20:200;
    
elseif row_num == 2 && col_num == 1
    d        = 8;
    N_max    = 2000;
    m_values = 50:50:500;
    
elseif row_num == 2 && col_num == 2
    d        = 16;
    N_max    = 4500;
    m_values = 100:100:1000;
    
else
    error('invalid subfigure selected');
end

num_m = length(m_values);

%%% Construcft index set and define N %%%

n = find_order(index_type,d,N_max);     % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n);       % compute index set
N = size(I,2);                          % number of basis functions

% arrays for storing the errors

L2_error_data   = zeros(num_m,num_trials,num_cases);
Linf_error_data = zeros(num_m,num_trials,num_cases);
residual        = zeros(1,num_cases);

%%% Construct error grid %%%

if ir_domain == 1                           % sampling points for the error grid
    err_samp_type = 'uniform';
else
    err_samp_type = sprintf('uniform_ir_%d',ir_domain);
end

M        = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(err_samp_type,d,M);


%%% Main loop %%%

for l = 1:num_cases
    
    % Legendre, uniform
    if l == 1
        if ir_domain == 1
            samp_type = 'uniform';
        else
            samp_type = sprintf('uniform_ir_%d',ir_domain);
        end
        basis_type     = 'legendre';        % functions used in the measurement matrix A
        err_basis_type = basis_type;        % polynomials used in computing the error
        
        % Legendre, optimal
    elseif l == 2
        samp_type      = 'optimal';
        basis_type     = 'legendre';
        err_basis_type = basis_type;
        
        % QR-Legendre, uniform
    elseif l == 3
        samp_type      = 'uniform_grid';
        basis_type     = 'orthogonal';
        err_basis_type = basis_type;
        I              = reorder_index_set(I,[]);          % Reorder index set
        
        % QR-Legendre, optimal
    else
        samp_type      = 'optimal';
        basis_type     = 'orthogonal';
        err_basis_type = basis_type;
        I              = reorder_index_set(I,[]);          % Reorder index set
    end
    
    %%% Construct error matrix and error vector %%%
    
    A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);
    
    if isequal(basis_type,'orthogonal')
        [Q_err_grid,R] = qr(A_err_grid,0);
        A_err_grid     = Q_err_grid;
    end
    
    b_err_grid = func(err_grid)/sqrt(M);
    
    %%% set the weights %%%
    
    if weighted_opt == 1        % weights via QR
        if isequal(basis_type,'orthogonal')
            w = max(abs(Q_err_grid))';
        else
            [Q_err_grid,R] = qr(A_err_grid,0);
            w = max(abs(Q_err_grid))';
        end
    else                        % non-weighted
        w = [];
    end
    
    % Compute optimal measure
    if isequal(samp_type,'optimal')
        [opt_mu,opt_w,theta] = opt_weights(A_err_grid*sqrt(M));
    else
        opt_mu = [];
        opt_w  = [];
    end
    
    %%% Compute a reference solution via least squares on the error grid %%%
    
    c_ref = A_err_grid\b_err_grid;
    residual(l) = norm(A_err_grid*c_ref-b_err_grid);
    
    for i = 1:num_m
        
        m = m_values(i);
        
        % temporary arrays to allow parfor
        L2_error_temp   = zeros(1,num_trials);
        Linf_error_temp = zeros(1,num_trials);
        
        %for t = 1:num_trials
        parfor t = 1:num_trials
            
            % generate sample points
            [y_grid, I_set] = generate_sampling_grid(samp_type,d,m,opt_mu,err_grid);
            
            % generate weights
            if isequal(samp_type,'optimal')
                W = diag(sqrt(opt_w(I_set)));
            else
                W = eye(m);
            end
            
            if isequal(samp_type,'uniform') || isequal(samp_type,'uniform_ir_2') || isequal(samp_type,'uniform_ir_3')
                % generate measurement matrix
                A = W*generate_measurement_matrix(basis_type,I,y_grid);
                
                % generate measurement vector
                b = W*func(y_grid)/sqrt(m);
            else
                % generate measurement matrix
                A = W*A_err_grid(I_set,:)*sqrt(M/m);
                
                % generate measurement vector
                b = W*b_err_grid(I_set)*sqrt(M/m);
            end
            
            % use reference solution to find best eta
            eta = norm(A*c_ref-b);
            
            % solve the weighted QCBP problem
            [c,stat] = wqcbp_cvx(A,b,w,eta,cvx_opt);
            
            % compute L^2_rho-norm and L^inf error
            L2_err   = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid);
            Linf_err = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf);
            
            disp(['Figure ',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),space,'dimension = ',num2str(d),space,...
                'polynomials = ',basis_type,space,'sampling = ',samp_type,space,'m = ',num2str(m),space,'trial = ',num2str(t),...
                space,'iter = ',num2str(i),space,'L2 error = ', num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
            
            L2_error_temp(1,t)   = L2_err;
            Linf_error_temp(1,t) = Linf_err;
            
        end
        
        L2_error_data(i,:,l)   = L2_error_temp;
        Linf_error_data(i,:,l) = Linf_error_temp;
        
    end
    
end

%%% Save data %%%
clear A_err_grid R b_err_grid err_grid Q_err_grid W w A b eta opt_mu opt_w theta

save(['../../data/Figure 11_12_13/',fig_name])


end