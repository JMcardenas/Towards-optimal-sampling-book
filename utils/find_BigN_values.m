%-------------------------------------------------------------------------%
%
% Filename: abs_rec.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description:  computes the reciprocal of summatory of absolute value and 
% square root of y
%
% Inputs:
% N_values - sizes of index set
% n_values - maximum degree of index set
% d - dimension
% opt - 1 input is N_values, 2 input is n_values
%
% Output:
% n_vals - n values, maximum degree of index sets 
% N_vals - N values, size of index sets
%-------------------------------------------------------------------------%

function [n_vals,N_vals] = find_BigN_values(N_values,n_values,d,opt)
 
if opt == 1
    % Input N_values (approximation of N_values)
    n_vals = zeros(1,length(N_values));
    N_vals = zeros(1,length(N_values));

    for i = 1:length(N_values)
        n_vals(1,i) = find_order('HC',d,N_values(i));    
        N_vals(1,i) = size(generate_index_set('HC',d,n_vals(1,i)),2);

    end
    
elseif opt == 2 
    % Input n_values
    N_vals = zeros(1,length(n_values));
    
    for i = 1:length(n_values)  
        N_vals(1,i) = size(generate_index_set('HC',d,n_values(i)),2);
        
    end
    n_vals = [];
    
else
   disp('wrong option N values') 
   
end
