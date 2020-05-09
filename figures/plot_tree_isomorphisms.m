function sp = plot_tree_isomorphisms(redundancy)

    fig_width  = 4.275;
    fig_height = 5;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    plot(redundancy, 'Color', 'blue', 'LineWidth', 1.0);
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    ylabel('Frequency', 'FontSize', 6);
        
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('Tree Rank', 'FontSize', 6);
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    title('Tree Isomorphisms', 'FontSize', 6, 'FontWeight', 'normal');
    
    axis tight;
    box off;
    
    left_margin = 0.6;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    right_margin = 0.4;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin]);
    
end