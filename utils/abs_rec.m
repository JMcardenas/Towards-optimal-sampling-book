%--- Description ---%
%
% Filename: abs_rec.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: computes the absolute reciprocal function f_4 
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = abs_rec(y)

%[m,d] = size(y); 
b     = 1./sum(sqrt(abs(y)),2);

end