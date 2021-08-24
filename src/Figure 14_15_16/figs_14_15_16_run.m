%-------------------------------------------------------------------------%
% Filename: figs_14_15_16_run.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates all data for Figures 14, 15 and 16
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath ../../utils

for fig_num = 14:16
    for row_num = 1:2
        for col_num = 1:2
            
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure ',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            figs_14_15_16_data(fig_num,row_num,col_num);
            
        end
    end
end