%--- Description ---%
%
% Filename: plot_book_style.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: Given random data associated 
% The plot shows the mean curve and a filled-in area defined by [mean - std
% dev, mean + std dev]. % For each curve, the data is assumed to be 
% relative to the same set of x values and to be generated using a fixed 
% number of random trials for each (x,y) value.
% 
% Inputs: 
% x_data - x-values (they should be the same for all curves)
% y_data - 3D array containing the y-values defiing the curves. This array
% should be structured as follows: y_data(i, j, k) contains the y-value 
% corresponding to the i-th x-value (i.e., x(i)), to the j-th random trial, 
% and the k-th curve.
% type - specifies what type of plot is produced.
%        The available options are:
%        'shaded'       - shaded plot (random data)
%        'shaded_pairs' - shaded plot (random data) with pairwise comparison 
%        'shaded_data'  - shaded plot (random data) + raw data
%        'line'         - line plot (deterministic data)
% stat - specifies what statistics of the data is visualized.
%        The available options are:
%        'mean_std_log10' - 10^( mean(log10(data)) +/- std(log10(data)) )
%        'mean_std'       - mean +/- std dev
%        'mean_std_eps'   - max(mean +/- std dev, machine precision)
%        'mean_sem'       - mean +/- standard error of the mean (SEM) 
% Note: If type corresponds to deterministic data, the value of stat is 
%       ignored.
%
% Output:
% hMeanPlots - handle correspoing to the mean plots only (to be used to add
% a legend)

function [hMeanPlots] = plot_book_style_each_curve(x_data, y_data, type, stat, n_curve)

%% Set data parameters
x_data = x_data(:)'; % make sure x is a column vector

% Extract y_data dimensions
n_x     = size(y_data,1);
n_trial = size(y_data,2);
%n_curve = size(y_data,3);
%n_curve = 1;

%% Set plotting parameters
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


%% Plot
switch type
    case'line'
        error('Not implemented')
    otherwise
        
        hMeanPlots = []; % initialize legend function handle
        for i_curve = n_curve : n_curve
            
            % Extract data realtive to current curve
            y_data_curve = y_data(:,:,1);
            
            % Define curves bounding the shaded region
            switch stat
                case 'mean_std_log10'
                    % compute statistics on log10 of data, i.e.
                    % 10^( mean(log10(data)) +/- std(log10(data)) )
                    mean_curve = 10.^(mean(log10(y_data_curve),2))';
                    std_curve = std(log10(y_data_curve),0,2)';
                    curve_min = 10.^(log10(mean_curve) - std_curve);
                    curve_max = 10.^(log10(mean_curve) + std_curve);
                case 'mean_std'
                    % mean +/- std dev
                    mean_curve = mean(y_data_curve,2)';
                    std_curve = std(y_data_curve,0,2)';
                    curve_min = mean_curve - std_curve;
                    curve_max = mean_curve + std_curve;
                case 'mean_std_eps'
                    % max(mean +/- std dev, machine precision)
                    mean_curve = mean(y_data_curve,2)';
                    std_curve = std(y_data_curve,0,2)';
                    curve_min = max(mean_curve - std_curve, eps);
                    curve_max = mean_curve + std_curve;
                case 'mean_sem'
                    % mean +/- standard error of the mean (SEM)
                    mean_curve = mean(y_data_curve,2)';
                    std_curve = std(y_data_curve,0,2)'/sqrt(n_trial);
                    curve_min = mean_curve - std_curve;
                    curve_max = mean_curve + std_curve;
                otherwise
                    error('Not implemented')
            end
            
            % Plot curves bounding the shaded region
            plot(x_data, curve_min, 'color', [colors{i_curve}, AlphaLevel]);
            hold on
            plot(x_data, curve_max, 'color', [colors{i_curve}, AlphaLevel])
            
            % Fill in shaded region
            x_data_2 = [x_data, fliplr(x_data)];
            inBetween = [curve_min, fliplr(curve_max)];
            fill(x_data_2, inBetween, colors{i_curve}, 'FaceAlpha', AlphaLevel, 'EdgeAlpha', 0);
            
            %dim_curve = [2 3 5 10 16];
            %display_name = sprintf('$d=%d$',dim_curve(n_curve));
            
            %display_name = sprintf('$d=2$');
            % Plot mean curve
            h = plot(x_data, mean_curve, markers{i_curve},...
                'markersize',ms,'MarkerFaceColor',colors{i_curve},...
                'MarkerEdgeColor',colors{i_curve},'LineWidth',lw,...
                'Color',colors{i_curve});
            %,'DisplayName',display_name);
            %legend(display_name)
            hMeanPlots = [hMeanPlots, h]; % add current mean plot to legend handle
        end 
        hold on
        
end
