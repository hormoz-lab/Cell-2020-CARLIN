function sp = plot_replicate_proportions(allele_breakdown_by_sample, is_template, title_str, show_legend)

    assert(size(allele_breakdown_by_sample,2)==2, 'Can only make replicate plot for two samples');
    
    if (nargin == 2)
        show_legend = false;
    end

    fig_width  = 4.4;
    fig_height = 3.0;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);       
    sp = subplot(1,1,1);
   
    denom = sum(allele_breakdown_by_sample,1);
    
    % Don't plot the template, because it is an outlier point, that crowds
    % the other points.
    
    allele_breakdown_by_sample = allele_breakdown_by_sample(~is_template,:);
    p = log10(allele_breakdown_by_sample./sum(allele_breakdown_by_sample,1));
    
    minval = min(p(isfinite(p)))-0.2;
    
    both = sum(logical(allele_breakdown_by_sample),2)==2;
    just1 = p(~both & allele_breakdown_by_sample(:,1)>0,1);
    just2 = p(~both & allele_breakdown_by_sample(:,2)>0,2);
    
    rng(23490);
    jitter_both = 0.1*rand(sum(both),2);
    jitter1 = 0.1*rand(size(just1));
    jitter2 = 0.1*rand(size(just2));    
    
    hold on;    
    scatter(p(both,1)+jitter_both(:,1), p(both,2)+jitter_both(:,2), 6, 'b', 'filled');
    scatter(just1+jitter1, minval*ones(length(just1),1), 12, 'r', 'filled');
    scatter(minval*ones(length(just2),1), just2+jitter2, 12, 'r', 'filled');
    hold off;
    
    axis tight;
    
    box off;
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'XAxis'), 'FontSize', 5);
        
    xlabel('log_{10} Allele Fraction in Rep 1', 'FontSize', 6);
    ylabel('log_{10} Allele Fraction in Rep 2', 'FontSize', 6);
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');

    title(title_str, 'interpreter', 'tex', 'FontSize', 6, 'FontWeight', 'normal');
    
    left_margin = 0.8;
    top_margin = 0.4;    
    bottom_margin = 0.7;        
    right_margin = 0.1;   
    
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin]);
        
    if (show_legend)
        lgd_height = 1.0;
        lgd_width = 1.8;
        inner_margin = 0.1;
        lgd = legend({sprintf('Both\nReplicates'); sprintf('Single\nReplicate')}, 'Location', 'Northwest', 'FontSize', 5);
        lgd.Units = 'centimeters';
        lgd.Position = [fig_width-inner_margin-right_margin-lgd_width bottom_margin lgd_width lgd_height];
        legend('boxoff');
    end   
end