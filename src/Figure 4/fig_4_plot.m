%-------------------------------------------------------------------------%
% Filename: figs_4_plot.m
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for learning sparse
% approximations in high dimensions", Springer, 2021.
%
% Description: generates the plot for Figure 4
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))


%% Set plotting parameters
type  = 'shaded';
stats = 'mean_std_log10';

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
%%

max_div     = 10;
num_trials  = 1;

Theta_new_data     = zeros(max_div,num_trials,5);
Big_Theta_new_data = zeros(max_div,num_trials,5);
N_values_data      = zeros(max_div,5);

for row_num = 1:2
    for col_num = 1:3
        
        if col_num == 1
            dim_num = 3;
        else
            dim_num = 5;
        end
        
        % load data
        for dim = 1:dim_num
            
            fig_load_name  = ['fig_4','_',num2str(row_num),'_',num2str(col_num),'_',num2str(dim)];
            load(['../../data/Figure 4/',fig_load_name])
            
            Theta_new_data(:,:,dim)     = Theta_data;
            Big_Theta_new_data(:,:,dim) = Big_Theta_data;
            N_values_data(:,dim)        = N_values;
        end
        
        if row_num == 1
            plot_data = Big_Theta_new_data;
        elseif row_num == 2
            plot_data = Theta_new_data;
        else
            disp(['error row num']);
        end
        
        fig = figure();
        
        for i_curve = 1:dim_num
            plot(N_values_data(:,i_curve),plot_data(:,:,i_curve),...
                markers{i_curve},'markersize',ms,'MarkerFaceColor',colors{i_curve},...
                'MarkerEdgeColor',colors{i_curve},'LineWidth',lw,'Color',colors{i_curve});
            hold on
        end
        
        if col_num == 1
            leg = legend('$d=1$','$d=2$','$d=4$','location','northwest');
        else
            leg = legend('$d=1$','$d=2$','$d=4$','$d=8$','$d=16$','location','northwest');
        end
        
        axis tight
        set_axis_param
        
        if row_num == 2
            ax = gca;
            ax.YTick = [0 5 10 15 20 25 30];
            
            if col_num == 1
                ylim([0 30]);
            elseif col_num == 3
                ylim([0 20])
            else
                disp('col num 1')
            end
        end
        
        set_fonts
        hold off
        fig_name_last  = ['fig_4','_',num2str(row_num),'_',num2str(col_num)];
        
        saveas(fig,['../../figs/Figure 4/',fig_name_last],'epsc');
        
    end
end
