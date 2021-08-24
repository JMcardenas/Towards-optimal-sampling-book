%-------------------------------------------------------------------------%
% Filename: fig_2_plot.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figure 2
%-------------------------------------------------------------------------%
% row_num - row number (either 1 or 2)
% col_num - column number (either 1, 2 or 3)
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all figures and subfigures %%%
for row_num = 1:2
    for col_num = 1:3
        
        % load name
        fig_name = ['fig_2','_',num2str(row_num),'_',num2str(col_num)];
        load(['../../data/Figure 2/',fig_name])
        
        fig = figure;
        
        plot(err_grid(:,1),err_grid(:,2),'.','Color',0.7*[1 1 1],'Markersize',4);
        hold on
        plot(err_grid(I_set,1),err_grid(I_set,2),'r.','Markersize',11);
        
        set_fonts
        set_axis_param
        ax = gca;
        grid off
        ax.XTick = [-1 -0.5 0 0.5 1];           xlim([-1 1]);
        ax.YTick = [-1 -0.5 0 0.5 1];           ylim([-1 1]);
        
        hold off
        
        fig_name_last = ['fig_2','_',num2str(row_num),'_',num2str(col_num)];
        saveas(fig,['../../figs/Figure 2/',fig_name_last],'epsc');
        
    end
end