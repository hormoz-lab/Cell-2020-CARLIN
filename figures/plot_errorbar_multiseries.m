function plot_errorbar_multiseries(dat, plt, wide)

    if (nargin == 2)
        wide = false;
    end

    N_series = unique(structfun(@(x) length(x), dat));
    assert(length(N_series)==1);

    fig_height = 5;
    if (wide)
        fig_width = 10;
    else
        fig_width  = 4.275;    
    end

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
        
    hold on;
    for i = 1:N_series
        if (isfield(dat, 'yerr'))
            errorbar(dat.xvals{i}, dat.yvals{i}, dat.yerr{i}, 'Marker', dat.markers{i}, 'DisplayName', dat.series_names{i}, ...
                 'Color', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 4, 'LineWidth', 1.0);
        else
            plot(dat.xvals{i}, dat.yvals{i}, 'Marker', dat.markers{i}, 'DisplayName', dat.series_names{i}, ...
                 'Color', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 4, 'LineWidth', 1.0);
        end
    end
    hold off;
    
    set(gca, 'LineWidth', 1.0);
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    
    if (isfield(plt, 'xticks'))
        xticks(plt.xticks); 
    end
    if (isfield(plt, 'xticklabels'))
        xticklabels(plt.xticklabels);
    end    
    if (isfield(plt, 'xlabel'))
        xlabel(plt.xlabel, 'FontSize', 6);
    end
        
    set(get(gca, 'YAxis'), 'FontSize', 5);
    if (isfield(plt, 'ylabel'))
        ylabel(plt.ylabel, 'FontSize', 6);
    end
    
    if (isfield(plt, 'xlim'))
        xlim(plt.xlim);
    end
    
    if (isfield(plt, 'ylim'))
        ylim(plt.ylim);
    end
    
    if (isfield(plt, 'title'))    
        title(plt.title, 'FontWeight', 'normal', 'FontSize', 6);
    end
    
    set(gca, 'TickDir', 'out');
    
    box off;
    
    legend('Location', 'SouthOutside', 'NumColumns', 4, 'FontSize', 6);
    legend('boxoff');

end