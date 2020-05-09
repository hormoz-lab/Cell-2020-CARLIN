function sp = plot_subtree_stability(instances, population, depth, N_sim)

    fig_width  = 5.7;
    fig_height = 5;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    line_color = lines;
        
    hold on;    
    plot(instances, 'color', line_color(1,:), 'LineWidth', 1.0);    
    ylim([0, N_sim*1.05]);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    ylabel('Frequency in Simulations', 'FontSize', 6);
    yyaxis right;
    plot(log10(population.median), 'color', line_color(2,:), 'LineWidth', 1.0);        
    fill([1:length(instances) length(instances):-1:1], ...
         log10([population.min; flipud(population.max)]), line_color(2,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');    
    yrulers = get(gca, 'YAxis');
    set(yrulers(1), 'Color', line_color(1,:));
    set(yrulers(2), 'FontSize', 5, 'Color', line_color(2,:));
    ylabel('log_{10} Clade Population', 'FontSize', 6);
    
    set(get(gca, 'XAxis'), 'FontSize', 5);        
    xlabel('Rooted Path Rank', 'FontSize', 6);
    
    hold off;
    axis tight;
    yl = ylim;
    ylim([0, yl(end)]);
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    title(sprintf('Stability of Depth-%d Rooted Paths', depth), 'FontSize', 6, 'FontWeight', 'normal');
    
    box off;
    
    left_margin = 0.9;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    right_margin = 0.8;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin]);
    
end