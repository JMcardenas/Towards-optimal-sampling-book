  
%-------------------------------------------------------------------------%
%
% Filename: F_Genz.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description:  computes the genz peak product function
%
% Inputs:
% y - input vector 
%
% Output:
% output - genz peak product function evaluated in y
%-------------------------------------------------------------------------%

function output = F_Genz(y)

[m,d] = size(y);
ind    = [1:d];
term_1 = (y+((-1).^(ind+1))./(ind+1)).^2;
term_2 = (d/4)./((d/4)+term_1);
output = prod(term_2,2);

end
