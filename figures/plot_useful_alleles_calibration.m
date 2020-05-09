function sp = plot_useful_alleles_calibration(N_cells, N_sig_alleles, N_clamp)

    assert(length(N_cells) == length(N_sig_alleles));

    fig_width  = 4.275;
    fig_height = 5;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    line_color = lines;
    
    N_peak = N_cells(N_sig_alleles == max(N_sig_alleles));
    max_cells = 1200;
    assert(max(N_sig_alleles) < max_cells);
    
    N_cells = [N_cells; N_clamp];
    N_sig_alleles = [N_sig_alleles; 0];
    N_sig_alleles(N_cells > N_clamp) = 0;
    
    [N_cells, reorder] = sort(N_cells);
    N_sig_alleles = N_sig_alleles(reorder);
        
    hold on;        
    patch([0 0 N_peak N_peak], [0 max_cells max_cells 0], 'green', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch([N_peak N_peak N_clamp N_clamp], [0 max_cells max_cells 0], 'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch([N_clamp N_clamp N_cells(end) N_cells(end)], [0 max_cells max_cells 0], 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(N_cells, N_sig_alleles, 'color', line_color(1,:), 'LineWidth', 1.0);        
    ylim([0, max_cells]);
    set(get(gca, 'YAxis'), 'FontSize', 5, 'Color', line_color(1,:));
    ylabel('Cells with Rare Alleles (\alpha=0.05)', 'FontSize', 6);
    yyaxis right;
    plot(N_cells, N_sig_alleles./N_cells, 'color', line_color(2,:), 'LineWidth', 1.0);    
    yrulers = get(gca, 'YAxis');    
    set(yrulers(1), 'Color', line_color(1,:));    
    set(yrulers(2), 'FontSize', 5, 'Color', line_color(2,:));    
    ylabel('Fraction of Total Cells', 'FontSize', 6);    
    
    set(get(gca, 'XAxis'), 'FontSize', 5);        
    xlabel('Edited Cells', 'FontSize', 6);
    axis tight;
    ylim([0 1]);    
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
  
    box off;
    
    left_margin = 0.8;
    top_margin = 0.1;    
    bottom_margin = 0.5;
    right_margin = 0.8;    
    leg_height = 1.0;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin+leg_height, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin-leg_height]);

    lgd = legend({'Optimal'; 'Saturated'; 'Undefined'}, ...
                 'Location', 'Northeast', 'FontSize', 5);
    lgd.Units = 'centimeters'; 
    lgd.NumColumns = 2;
    
    lgd.Position = [0, 0, fig_width, leg_height];
    
    legend('boxoff');

end