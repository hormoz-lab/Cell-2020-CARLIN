function sp = plot_allele_on_transcriptome(experiment, x, y, ceres, ...
                                           ceres_names, sz, color, marker, plot_order_bg_to_fg, ...
                                           show_legend, legend_order, title_str, mini, highlight_patch)
                                           
    assert(length(x) == length(ceres) && length(y) == length(ceres));
    N_series = size(ceres_names,1);
    assert(max(ceres) <= N_series);
    assert(length(sz) == N_series);
    assert(size(color,1) == N_series);
    assert(isequal(sort(plot_order_bg_to_fg), [1:N_series]'));
    assert(isequal(sort(legend_order), [1:N_series]'));

    if (mini)
        fig_width  = 2.75;
        fig_height = 2.5;
    else
        fig_width  = 4.7;
        fig_height = 4.3;
    end
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    if(isscalar(sz))
        sz = repelem(sz, size(color,1));
    end

    hold on;       
    scatter(x, y, 3, [0.9 0.9 0.9], 'filled');     
    if (~isempty(highlight_patch))
        fill(highlight_patch(:,1), highlight_patch(:,2), [0 0.7 0], 'EdgeColor', [0 0.7 0], ...
             'FaceAlpha', 0.2, 'EdgeAlpha', 0.5, 'LineWidth', 1.0, 'DisplayName', 'HSC Cluster');
    end
    
    [~, reorder] = sort(plot_order_bg_to_fg);
    
    for i = reorder'
        if (any(ceres==i))            
            scatter(x(ceres==i), y(ceres==i), sz(i), color(i,:), 'filled', marker(i), 'DisplayName', ceres_names{i});            
        else
            scatter(NaN, NaN, sz(i), color(i,:), 'filled', marker(i), 'DisplayName', ceres_names{i});
        end
    end
    hold off;
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal', 'interpreter', 'none');
    dat_series = get(gca, 'Children');
    scatter_series = flipud(dat_series(1:length(reorder)));
    
    [~, reorder] = sort(reorder);
    scatter_series = scatter_series(reorder);
    scatter_series = scatter_series(legend_order);
    if (highlight_patch)
        dat_series = [dat_series(end-1); scatter_series];
    else
        dat_series = scatter_series;
    end
    
    axis tight;
    axis off;
    box off;

    tight_margin = 0.3;

    left_margin = 0.1;
    right_margin = 0.1;
    top_margin = 0.25;
    bottom_margin = 0.1;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-right_margin-left_margin fig_height-top_margin-bottom_margin]);
  
    if (show_legend)
        lgd = legend(dat_series, 'FontSize', 5);
        lgd.Units = 'centimeters';
        
        if (strcmp(experiment, 'EB'))
            leg_width = 2.5;
            leg_height = 0.6;
            lgd.Position = [tight_margin, fig_height-top_margin-tight_margin-leg_height, leg_width, leg_height];

        elseif (strcmp(experiment, '5FU'))
            leg_width = 2.5;
            leg_height = 1.2;
            lgd.Position = [tight_margin tight_margin, leg_width, leg_height];        
        end  
    end
end