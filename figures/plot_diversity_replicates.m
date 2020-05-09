function fig = plot_diversity_replicates(N_sampled, N_unique1, N_unique2, N_shared, N_union, pval)

    assert(size(N_unique1,1) == length(N_sampled));
    assert(size(N_unique2,1) == length(N_sampled));
    assert(size(N_shared, 1) == length(N_sampled));
    assert(size(pval, 1)     == length(N_sampled));
 
    assert(isequal(cellfun(@length, N_union), cellfun(@length, pval)));
    assert(isequal(cellfun(@length, N_unique1)+cellfun(@length, N_unique2)+cellfun(@length, N_shared), cellfun(@length, N_union)));
 
    fig_width  = 10.5;
    fig_height = 8;
    fig = figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
                 'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
             
    sp{1} = subplot(3,1,1);
    b = bar(1:length(N_sampled), [cellfun(@length, N_unique1(:,1)), cellfun(@length, N_unique2(:,1)), cellfun(@length, N_shared(:,1))], ...
            'stacked', 'BarWidth', 1.0);    
    b(1).FaceColor = [0 0 1];
    b(2).FaceColor = [0.5 0.5 1];
    b(3).FaceColor = 'r';
    axis tight;
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    ylabel('# Alleles', 'FontSize', 6);
    
    set(gca, 'xtick', []);
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    box off;
    
    lgd{1} = legend({'Mouse 1 Only'; 'Mouse 2 Only'; 'Both Mice'}, 'Location', 'Northwest', 'FontSize', 5);
    title(lgd{1}, 'Allele Membership', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
            
    sig_level = 0.05;
    horz_offset = 0.2;
    jitter_val = 0.1;
    
    sp{2} = subplot(3,1,2);
    hold on;
    
    h(1) = line([1-horz_offset-jitter_val length(N_sampled)+horz_offset+jitter_val], [sig_level, sig_level], 'LineStyle', '--', 'Color', 'black');
    tempx = repelem(1:length(N_sampled), cellfun(@(x,y) sum(~ismember(x,y)), N_union(:,1), N_shared(:,1)));
    tempy = cellfun(@(x,y,z) x(~ismember(y,z)), pval(:,1), N_union(:,1), N_shared(:,1), 'un', false);
    h(2) = scatter(tempx-horz_offset, vertcat(tempy{:}), 1, 'b', 'filled', 'jitter', jitter_val);
    tempx = repelem(1:length(N_sampled), cellfun(@(x,y) sum(ismember(x,y)), N_union(:,1), N_shared(:,1)));
    tempy = cellfun(@(x,y,z) x(ismember(y,z)), pval(:,1), N_union(:,1), N_shared(:,1), 'un', false);
    h(3) = scatter(tempx+horz_offset, vertcat(tempy{:}), 1, 'r', 'filled', 'jitter', jitter_val);
    hold off;
    axis tight;
    
    set(gca, 'yscale', 'log');
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    ylabel('p_{clonal}', 'FontSize', 6);
    
    set(gca, 'xtick', []);
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    box off;
    
    lgd{2} = legend(h([2 3]), {'1 Mouse'; 'Both Mice'}, 'Location', 'Southeast', 'FontSize', 5);    
    legend('boxoff');
        
    sp{3} = subplot(3,1,3);    
    frac_shared = cellfun(@(p,u,s) sum(p(ismember(u, s)) < sig_level)/length(s), pval, N_union, N_shared);
    frac_sig    = cellfun(@(p,u,s) sum(p(ismember(u, s)) < sig_level)/sum(p < sig_level), pval, N_union, N_shared);
    
    hold on;    
    errorbar(1:length(N_sampled), nanmean(frac_sig,2), nanstd(frac_sig, [], 2), 'b', 'LineWidth', 1.0);    
    errorbar(1:length(N_sampled), nanmean(frac_shared,2), nanstd(frac_shared, [], 2),  'r', 'LineWidth', 1.0);    
    line([1-horz_offset-jitter_val length(N_sampled)+horz_offset+jitter_val], [sig_level, sig_level], 'LineStyle', '--', 'Color', 'black');
    hold off;
    axis tight;
    ylim([0 0.2]);
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    ylabel('Fraction', 'FontSize', 6);
    
    set(gca, 'xtick', [1:length(N_sampled)]);
    set(gca, 'xticklabels', arrayfun(@num2str, N_sampled, 'UniformOutput', false));
    xtickangle(90);
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('Total Edited Cells Sampled', 'FontSize', 6);
    
    box off;
        
    lgd{3} = legend({'normalized by # of total alleles with p_{clonal} < 0.05'; 
                     'normalized by # of shared alleles'}, ...
                     'Location', 'Northeast', 'FontSize', 5);
    title(lgd{3}, '# of shared alleles with p_{clonal} < 0.05 ...', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');

    linkaxes([sp{1}; sp{2}; sp{3}], 'x');
    
    left_margin  = 1.1;
    right_margin = 0.1;
    bottom_margin = 1.0;
    top_margin = 0.1;
    inner_margin = 0.2;
    
    panel_height = (fig_height-bottom_margin-top_margin-2*inner_margin)/3;
    panel_width = fig_width-left_margin-right_margin;
    
    set(sp{1}, 'Units', 'centimeters', 'Position', [left_margin,  bottom_margin+2*panel_height+2*inner_margin, panel_width, panel_height]);
    set(sp{2}, 'Units', 'centimeters', 'Position', [left_margin,  bottom_margin+1*panel_height+1*inner_margin, panel_width, panel_height]);
    set(sp{3}, 'Units', 'centimeters', 'Position', [left_margin,  bottom_margin,                               panel_width, panel_height]);
    
    leg_width = 5.5;
    leg_height = 1;
    lgd{3}.Units = 'centimeters';        
    lgd{3}.Position = [fig_width-leg_width bottom_margin+panel_height-leg_height leg_width leg_height];
    
end