%-------------------------------------------------------------------------%
%
% Filename: const_theta.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for learning sparse
% approximation in high dimensions", Springer, 2021.
%
% Description:  computes the constant theta for measurement matrix A
%
% Inputs: 
% A_err - error matrix dimensions [size(err_grid), size(Index_set)]
%
% Output:
% theta - constant theta 
%-------------------------------------------------------------------------%

function theta = const_theta(A_err)

sup_psi = max((abs(A_err).^2),[],2);        
theta   = sqrt(mean(sup_psi));              % Integration via MC

end
