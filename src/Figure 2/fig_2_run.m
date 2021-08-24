%-------------------------------------------------------------------------%
% Filename: fig_2_run.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates all data for Figure 2
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath ../../utils

for row_num = 1:2
    for col_num = 1:3
        disp(' ');
        disp(['------------------------------------------------------------------------']);
        disp(['Running Figure 2','_',num2str(row_num),'_',num2str(col_num)]);
        disp(['------------------------------------------------------------------------']);
        disp(' ');
        fig_2_data(row_num,col_num);
    end
end
