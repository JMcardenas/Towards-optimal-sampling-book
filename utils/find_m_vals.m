%-------------------------------------------------------------------------%
%
% Filename: find_m_vals.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description:  computes the reciprocal of summatory of absolute value and 
% square root of y
%
% Inputs:
% d- dimension
% n_max - maximum degree of index set
% max_div - maximum division
% N_values - sizes of index sets 
% opt - option, 1 input is n_max, 2 input is N_values
%
% Output:
% m_values - m samples values 
% n_values - n values, degree of index sets
% N_values - size of index sets
% k_ratios - k ratios 
%-------------------------------------------------------------------------%

function [m_values,n_values,N_values,k_ratios] = find_m_vals(d,n_max,max_div,N_values,opt)


index_type = 'HC';

if opt == 1
    % Input is n_max
    N_values = zeros(1,max_div-1);
    n_values = linspace(1,n_max,max_div); 

    for j = 1:(max_div-1) 
        I = generate_index_set(index_type,d,n_values(j));   % compute index set
        N = size(I,2);

        N_values(1,j) = N;
    end

    k_ratios = round(log(N_values));
    m_values = k_ratios.*N_values;
    
elseif opt == 2
    % Input is N_values
    k_ratios = round(log(N_values));
    m_values = k_ratios.*N_values;
    n_values = [];
    
elseif opt == 3
    % Input is m_values
    
else
    disp('wrong option');
    
end


end
