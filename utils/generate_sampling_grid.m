%-------------------------------------------------------------------------%

% Filename: generate_sampling_grid.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description: generates a random sampling grid using either the uniform, 
% Chebyshev (arcsine) measure, or optimal measure defined by 7.5.3.
%
% Inputs: 
% samp_type - either 'uniform' (uniform measure) or 'Chebyshev' 
%    (Chebyshev measure)
% d - dimension
% m - number of sample points
%
% Output:
% y_grid - the m x d array where the ith row is the ith sample point y_i
%-------------------------------------------------------------------------%
 
function [y_grid, I_set] = generate_sampling_grid(samp_type,d,m,opt_mu,err_grid)

% 1st uniform measure
if isequal(samp_type,'uniform')
    
    y_grid = 2*rand(m,d)-ones(m,d); 
    I_set  = [];
    
% 2nd uniform measure
elseif isequal(samp_type,'uniform_grid')
    
    K      = size(err_grid,1);
    I_set  = datasample([1:K]',m,'Replace',true); 
    y_grid = err_grid(I_set,:);
    
% Chebyshev measure
elseif isequal(samp_type,'chebyshev')

    y_grid = rand(m,d);
    y_grid = cos(pi .* y_grid); 
     
% Optimal measure or Adaptive measure
elseif isequal(samp_type,'optimal') || isequal(samp_type,'adaptive') || isequal(samp_type,'nonadaptive')
    
    K      = size(err_grid,1);
    I_set  = datasample([1:K]',m,'Replace',true,'Weights',opt_mu); 
    y_grid = err_grid(I_set,:);
    
% uniform measure over Annulus
elseif isequal(samp_type,'uniform_ir_2')
    
    y_grid = [];
    I_set  = [];
    rad    = 1;
    %rad   = 1/2;
    
    while length(y_grid) < m
    
        u = randn(1,d);
        x = rad*(u/norm(u))*rand^(1/d);
        
        if norm(x)>rad/4
            y_grid = [y_grid ; x];
        end
    end
    
% uniform measure over Linear constraint
elseif isequal(samp_type,'uniform_ir_3')
    
    y_grid = [];
    I_set  = [];
    
    while length(y_grid) < m        
        
        x = 2*rand(1,d)-1;
        
        if sum(x) <= 1  
            y_grid = [y_grid ; x];
        end
    end
    
% uniform measure over Prism
elseif isequal(samp_type,'uniform_ir_4')
    
    y_grid = [];
    I_set  = [];
    
    while length(y_grid) < m
        
        x = 2*rand(1,d)-1;
        
        if norm([x(1),x(2)]) >= (1/2)
            y_grid = [y_grid ; x];
        end
    end
    
% Invalid samp_type    
else

error('invalid samp_type')

end
