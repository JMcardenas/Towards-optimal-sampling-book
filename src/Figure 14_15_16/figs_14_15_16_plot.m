%-------------------------------------------------------------------------%
% Filename: figs_14_15_16_data.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the data for Figures 14, 15 and 16
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all figures and subfigures %%%

for fig_num = 14:16
    for row_num = 1:2
        for col_num = 1:2
            
            fig_name = ['fig','_',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            load(['../../data/Figure 14_15_16/',fig_name])
            
            fig = figure;
            
            hPlot = plot_book_style(m_values, Linf_error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight
            
            legend(hPlot, 'QU','QO','QUw','QOw','location','northeast');
            
            set_axis_param
            set_fonts
            saveas(fig,['../../figs/Figure 14_15_16/',fig_name],'epsc');
        end
    end
end
