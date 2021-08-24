%-------------------------------------------------------------------------%
% Filename: figs_8_9_10_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates all data for Figures 8, 9 and 10
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath ../../utils

for fig_num = 8:10
    for row_num = 1:2
        for col_num = 1:3
            
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure ',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            figs_8_9_10_data(fig_num,row_num,col_num);
            
        end
    end
end
