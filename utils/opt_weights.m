%-------------------------------------------------------------------------%
%
% Filename: opt_weights.m
% Programmer: Juan M. Cardenas and Ben Adcock
% Part of the chapter book "Towards optimal sampling for sparse 
% approximation in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figures 7.5 and 7.6
%
% Inputs:
% A_err - error matrix dimensions [size(err_grid), size(Index_set)]
%
% Output: 
% opt_mu - optimal sampling measure
% opt_w  - optimal weights
% theta  - theta value 
%-------------------------------------------------------------------------%

 
function[opt_mu,opt_w,theta] = opt_weights(A_err)

% sup over the basis 
sup_psi = max((abs(A_err).^2),[],2);        

% theta
theta   = sqrt(mean(sup_psi));              % Integration via MC

% optimal measure
inv_w   = (theta^(-2))*sup_psi;
opt_mu  = inv_w/size(A_err,1);

% weights
opt_w = 1./inv_w;

end
