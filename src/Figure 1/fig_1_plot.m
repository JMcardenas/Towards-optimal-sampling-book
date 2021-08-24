%-------------------------------------------------------------------------%
% Filename: fig_1_plot.m
% Programmer: Juan M. Cardenas and Ben Adcock
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates the plots for Figure 1
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Set plotting parameters %%%

type = 'shaded';

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

% set style based on plot type
switch type
    case 'shaded'
        % colors do not need to be changed in this case
    case 'shaded_pairs'
        % adapt colors and markers for comparison plot
        for i = 1:length(colors)
            new_colors{2*i-1} = colors{i};
            new_colors{2*i}   = 0.7 * colors{i};
        end
        new_markers = markers;
        for i = 2:2:length(markers)
            new_markers{i} = ['-',markers{i}];
        end
        colors = new_colors;
        markers = new_markers;
    otherwise
        error('Not implemented')
end

%%% Plot %%%

% define parameters
num_m_approx = 11;
num_trials   = 50;
num_dim      = 5;

% auxiliary variables
L2_error_plot_data = zeros(num_m_approx,num_trials,num_dim);
m_values_plot_data = zeros(num_m_approx,num_dim);

%%% Loop over all figures and subfigures %%%
for row_num = 1:3
    for col_num = 1:3
        
        h_final = [];
        fig = figure;
        
        for dim= 1:num_dim
            
            % load data
            fig_name = ['fig_1','_',num2str(row_num),'_',num2str(col_num),'_',num2str(dim)];
            load(['../../data/Figure 1/',fig_name])
            
            % save data
            L2_error_plot_data   = L2_error_data;
            Linf_error_plot_data = Linf_error_data;
            
            % generate plot
            hPlot   = plot_book_style_each_curve(m_values, Linf_error_plot_data, 'shaded', 'mean_std_log10',dim);
            h_final = [h_final, hPlot];
            
            set(gca, 'yscale', 'log','xscale','log')
            axis tight
            
        end
        
        legend(h_final,'$d=2$','$d=3$','$d=5$','$d=10$','$d=16$','interpreter','latex','location','southwest');
        
        set_axis_param
        set_fonts
        
        ax = gca;
        ax.XTick = [1e0 1e1 1e2 1e3 1e4];    xlim([1e0  1e4]);
        
        fig_name_last = ['fig_1','_',num2str(row_num),'_',num2str(col_num)];
        saveas(fig,['../../figs/Figure 1/',fig_name_last],'epsc');
        
    end
end






