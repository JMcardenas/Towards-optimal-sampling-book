%-------------------------------------------------------------------------%
% Filename: fig_3_plot.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates and saves the plots for Figure 3
%-------------------------------------------------------------------------%
% row_num - row number (either 1 or 2)
% col_num - column number (either 1, 2 or 3)
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all figures and subfigures %%%

for row_num = 1:2
    for col_num = 1:3
        
        load_name = ['fig_3','_',num2str(row_num),'_',num2str(col_num)];
        load(['../../data/Figure 3/',load_name])
        
        C_data_plot = zeros(num_m,num_trials,5);
        
        for l = 1:num_m_type
            C_data_plot(:,:,l) = C_data(:,:,1,l);
        end
        
        fig_name = ['fig_3','_',num2str(row_num),'_',num2str(col_num)];
        
        fig = figure();
        hPlot = plot_book_style(N_values, C_data_plot , 'shaded', 'mean_std_log10');
        grid on
        
        set(gca, 'yscale', 'log')
        axis tight
        
        legend(hPlot, '$m=s$','$m=2s$','$m=3s$','$m=s\log(s)$','$m=2s\log(s)$','location','northwest');
        
        ax = gca;
        ax.XTick = [0 200 400 600 800 1000];    xlim([0 1000]);
        ax.YMinorTick = 'off';          ax.YMinorGrid = 'off';
        ax.YGrid      = 'on';           ax.XGrid      = 'on';
        ax.XMinorTick = 'off';          ax.XMinorGrid = 'off';
        ax.YMinorGrid = 'off';          ax.YMinorTick = 'off';
        
        if col_num == 1
            ax.YTick = [1e0 1e5 1e10 1e15 1e20];
            ylim([1e0 1e17]);
            
        else
            ax.YTick = [1e0 1e1 1e2 1e3 1e4];
            ylim([1e0 1e3]);
        end
        
        %set_axis_param
        set_fonts
        saveas(fig,['../../figs/Figure 3/',fig_name],'epsc');
    end
end
