%-------------------------------------------------------------------------%
% Filename: figs_11_12_13_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates all data for Figure 11, 12 and 13
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all figures and subfigures %%%

for fig_num = 11:13
    for row_num = 1:2
        for col_num = 1:2
            
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure ',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            figs_11_12_13_data(fig_num,row_num,col_num);
            
        end
    end
end