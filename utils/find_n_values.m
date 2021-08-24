%-------------------------------------------------------------------------%
%
% Filename: find_n_vals.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description:  computes the reciprocal of summatory of absolute value and 
% square root of y
%
% Inputs:
% d- dimension
% N_max - maximum size of index set 
%
% Output: 
% n_values - n values, degree of index sets 
%-------------------------------------------------------------------------%

function n_values = find_n_values(N_max,d)
     
n_max = find_order('HC',d,N_max);

if d < 10
    switch d 
        case 2 
            r = 13;
        case 3 
            r = 13;
        case 5 
            r = 15;
    end
    
    n_val = 1:n_max;
    
    for i = 1 : length(n_val)
        In_t     = generate_index_set('HC',d,n_val(i));
        N_val(i) = length(In_t); 
        M_val(i) = round(log(N_val(i)))*N_val(i);
    end
    
    p1    = log(M_val(1));         
    p2    = log(M_val(end));
    M_p   = linspace(p1,p2,r);    
    pts   = exp(M_p); 
    j     = 1; 
    Index = []; 
        
    %--- Find the M_val closest to pts ---% 
    for i = 1 : r 
        while (M_val(j)<= round(pts(i)))  
            j = j + 1;
            if j == n_max + 1
                break  
            end
        end
        
        Index = [Index; j-1];
    end
    
    Index    = unique(Index);
    n_values = n_val(Index); 
    
else 
    n_values = 1:n_max;

end  


end
