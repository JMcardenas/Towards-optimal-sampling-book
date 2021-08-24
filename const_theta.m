%-------------------------------------------------------------------------%
%
% Filename: opt_weights.m
% Programmer: Juan M. Cardenas
% Part of the chapter book "Towards optimal sampling for sparse 
% approximation in high dimensions", Springer, 2021.
%
% Description: This function computes the optimal weights from sec 7.5.3
% 
% Inputs: 
% A_err - error matrix dimensions [size(err_grid), size(Index_set)]

function[theta] = const_theta(A_err)

sup_psi = max((abs(A_err).^2),[],2);        
theta  = sqrt(mean(sup_psi));              % Integration via MC

end
