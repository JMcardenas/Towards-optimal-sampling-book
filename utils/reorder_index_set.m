%-------------------------------------------------------------------------%
%
% Filename: reorder_index_set.m 
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description:   This function reorder the index set for the matrix Q.
%
% Inputs: 
% I - index set
%
% Output:
% J - new index set 
%-------------------------------------------------------------------------%

function J = reorder_index_set(I,n_values)

k     = max(sum(I));
[d,N] = size(I); 
J     = [];

if isempty(n_values) == 1                   % length(n_values) == 0
    % Re-order index set for QR using CS 
    for j = 0:k
        for i = 1:N
            if sum(I(:,i)) == j
                J = [J, I(:,i)];
            end
        end
    end
    
elseif length(n_values) > 1
    % Re-order index set for ASGD using LS 
    length_n = length(n_values); 
   
    for i = 1 : length_n

        I  = generate_index_set('HC',d,n_values(i));

        if i == 1
            J = I;
        else
            C = setdiff(I',J','rows');              % Index in I\J
            C = C';                                   
            J = [J  C];                             % adding K to J indices
        end
    end 
else
    disp('incorrect reorder index set ')
end
