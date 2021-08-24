%-------------------------------------------------------------------------%
% Filename: figs_8_9_10_data.m 
% Authors: Juan M. Cardenas and Ben Adcock.
% Part of the book chapter "Towards optimal sampling for sparse 
% approximation in high dimensions", Springer, 2021. 
%
% Description: generates the plots for Figures 8, 9 and 10
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../../utils'))
addpath(genpath('../../../graphics'))
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%%% Loop over all figures and subfigures %%% 
for fig_num = 8:10
    for row_num = 1:2
        for col_num = 1:3       

            fig_name = ['fig','_',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            load(['../../data/Figure 8_9_10/',fig_name])

            fig = figure();

            for j=1:2
                semilogy(coeff_sort(:,j),'o','markersize',ms,'MarkerFaceColor',...
                    colors{j},'MarkerEdgeColor',colors{j},'LineWidth',lw,...
                    'Color',colors{j});
                hold on
            end
            hold off 

            axis tight
            legend('$\Phi$','$\Upsilon$','interpreter','latex','location','northeast'); 

            set_axis_param
            set_fonts
            ax = gca;
            ax.YTick = [1e-15 1e-10 1e-5 1e0];
            ylim([1e-15 1e0])


            saveas(fig,['../../figs/Figure 8_9_10/',fig_name],'epsc');
            
        end
    end
end
