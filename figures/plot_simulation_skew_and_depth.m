function plot_simulation_skew_and_depth(prolif_b, sample_fact, dup_prob_induced, dup_prob_sampled)
        
    [N_beta, N_f] = size(dup_prob_sampled);
    assert(N_beta == length(prolif_b) && N_f == length(sample_fact));

    fig_width  = 7.0;
    fig_height = 3.9;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
     
    p = plot(log10(prolif_b), dup_prob_sampled/dup_prob_induced, 'marker', 'o', 'MarkerSize', 3, 'LineWidth', 1.0);
    set(p, {'MarkerFaceColor'}, get(p,'Color'));
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    xlabel('log_{10} Proliferation Probabilility Skew (b)', 'FontSize', 6);
    ylabel('$\gamma(f)$', 'FontSize', 6, 'interpreter', 'latex');
    axis tight;
    ylim([0 1.1]);
    title('Fraction of Multi-mapped Alleles with Multiple Observed Progenitors', 'FontSize', 6, 'FontWeight', 'normal');
    lgd = legend(cellstr(num2str(sample_fact)), 'FontSize', 5, 'Location', 'Northeast');
    title(lgd, 'Sampling Depth, f', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    box off;
    
    left_margin = 0.7;
    top_margin = 0.3;    
    bottom_margin = 0.8;
    right_margin = 0.1;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin]);
    
end

    