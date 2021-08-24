%-------------------------------------------------------------------------%
% Filename: fig_1_run.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates all data for Figure 1.
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath ../../utils

for row_num = 1:3
    for col_num = 1:3
        for dim = 1:5
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure 1','_',num2str(row_num),'_',num2str(col_num),'_',num2str(dim)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            fig_1_data(row_num,col_num,dim);
        end
    end
end

