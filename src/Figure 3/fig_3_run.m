%-------------------------------------------------------------------------%
% Filename: fig_3_run.m 
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse 
% approximations in high dimensions", Springer, 2021. 
%
% Description: generates and saves the data for Figure 3
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath ../../utils
 
for row_num = 1:2
    for col_num = 1:3
        disp(' ');
        disp(['------------------------------------------------------------------------']);
        disp(['Running Figure 3','_',num2str(row_num),'_',num2str(col_num)]);
        disp(['------------------------------------------------------------------------']);
        disp(' ');
        fig_3_data(row_num,col_num);  
    end
end

