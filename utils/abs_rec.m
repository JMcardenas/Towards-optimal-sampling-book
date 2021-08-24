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
% y - input vector 
%
% Output:
% b - function evaluated in y
%-------------------------------------------------------------------------%

function b = abs_rec(y)

b = 1./sum(sqrt(abs(y)),2);

end
