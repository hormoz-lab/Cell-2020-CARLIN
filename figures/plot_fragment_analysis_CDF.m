function sp = plot_fragment_analysis_CDF(frag, algo, L, envelope, title_str, show_legend)

    fig_width  = 5.8;
    fig_height = 2.5;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    hold on;
    plot(frag.xvals, frag.CDF, 'b', 'LineWidth', 1.0);
    plot(algo.xvals, algo.CDF, 'r', 'LineWidth', 1.0);
    patch([L-envelope(1) L-envelope(1) L+envelope(2) L+envelope(2)], [0 1 1 0], ...
          'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none');    
    hold off;
    axis tight;
        
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    ylabel('CDF', 'FontSize', 6);
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('Allele Length (bp)', 'FontSize', 6);
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    box off;
    
    left_margin = 0.7;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    tight_margin = 0.2;        
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-tight_margin fig_height-top_margin-bottom_margin]);
    
    if (show_legend)
        lgd = legend({sprintf('Fragment Analyzer'); sprintf('CARLIN Pipeline'); sprintf('Unedited Allele')}, ...
                     'Location', 'Northwest', 'FontSize', 5);
        lgd.Units = 'centimeters';        
        leg_width = 2.0;
        leg_height = 1.0;
        lgd.Position = [left_margin+3*tight_margin fig_height-top_margin-tight_margin-leg_height, leg_width, leg_height];
        legend('boxoff');
    end
    
end