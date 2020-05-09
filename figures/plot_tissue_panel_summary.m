function plot_tissue_panel_summary(tissue_labels, yvals, title_str, yl, series_names, body)

    if (body)        
        fig_width  = 4.275;
        fig_height = 5;        
    else        
        N_cells = cellfun(@length, tissue_labels);
        tissue_labels = vertcat(tissue_labels{:});
        yvals = vertcat(yvals{:});
         
        end_cells = cumsum(N_cells);
        start_cells = [1; end_cells(1:end-1)+1];
        end_cells   = end_cells+0.4;
        start_cells = start_cells-0.4;
        
        fig_width  = 6.8;
        fig_height = 5;
    end

    [N_tissues, N_samples] = size(yvals);
    assert(N_samples <= 3);
    
    assert(length(tissue_labels) == N_tissues);
    assert(length(series_names)  == N_samples);
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);

    left_margin = 0.7;
    top_margin = 0.3;    
    bottom_margin = 1.0;
    right_margin = 0.1;
    
    colors = ['rb'];
    markers = ['^'; 's'];
        
    hold on;    
    for i = 1:N_samples
        scatter([1:N_tissues], yvals(:,i), 12, colors(i), 'filled', 'Marker', markers(i));
    end
    
    if (~body)
        for i = 1:length(N_cells)
            line([start_cells(i) end_cells(i)], [100 100], 'Color', 'black', 'LineWidth', 2.0);
        end
    end
    
    hold off;
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(gca, 'Xtick', [1:N_tissues], 'xticklabel', tissue_labels);
    xtickangle(45);
    
    axis tight;
    ylim([0, 100]);
    
    ylabel(yl, 'FontSize', 6);
    
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    if (body)
        lgd = legend(series_names, 'FontSize', 5, 'Location', 'Northwest');    
    else
        lgd = legend(series_names, 'FontSize', 5, 'Location', 'Northeast');    
    end
    legend('boxoff');
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
end 