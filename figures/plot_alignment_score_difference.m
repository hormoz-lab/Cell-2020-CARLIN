function sp = plot_alignment_score_difference(score_diff, score_diff_CDF)

    assert(nargin == 2)
    
    fig_width  = 5.8;
    fig_height = 3.0;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);       
    sp = subplot(1,1,1);
    
    hold on;
    
    assert(all(score_diff >= 0));
    assert(min(score_diff_CDF) >= 0 && max(score_diff_CDF) <= 1);
    
    plot(score_diff, score_diff_CDF, 'Color', [0 0 1], 'LineWidth', 1.0);  
    axis tight;
    ylim([0 1]);
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('\Delta NW Score (NW Aligned - CARLIN Aligned)', 'FontSize', 6);
    ylabel('CDF (over alleles)', 'FontSize', 6);    
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    title('Relative NW Alignment Score', 'FontSize', 6, 'FontWeight', 'normal');
    
    box off;
    
    left_margin = 0.7;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    tight_margin = 0.1;        
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, ...
                                                 fig_width-left_margin-tight_margin fig_height-top_margin-bottom_margin]);
    
end