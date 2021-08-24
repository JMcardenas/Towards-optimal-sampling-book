%-------------------------------------------------------------------------%
% Filename: fig_4_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 4
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath ../../utils

for row_num = 1:2
    for col_num = 1:3
        
        if col_num == 1
            dim_num = 3;
        else
            dim_num = 5;
        end
        
        for dim = 1:dim_num
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure 4','_',num2str(row_num),'_',num2str(col_num),'_',num2str(dim)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            fig_4_data(row_num,col_num,dim);
        end
    end
end
