function make_invivo_phylogeny_subplots(results_dir)

    load(sprintf('%s/Trees/InVivoPhylogeny/Lineage.mat', results_dir));
    
    outdir = sprintf('%s/Figures/InVivoPhylogeny', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end 
    
    [pruned_adj_list, pruned_allele_at_node, N_pruned_alleles] ...
            = generate_adj_list_from_subtrees(stable_depthX_subtrees, depthX_subtree_clade_pop, 'min', 100);

    pruned_pop = sum(allele_breakdown_by_sample(pruned_allele_at_node(1:N_pruned_alleles),:),2);
    tissue_frac = allele_breakdown_by_sample./sum(allele_breakdown_by_sample,2);

    [~, reorder_hook] = arrayfun(@(i) sort(tissue_frac(pruned_allele_at_node(i),:), 'descend'), [2:N_pruned_alleles]', 'un', false);
    [~, reorder_hook] = sortrows(vertcat(reorder_hook{:}));
    [~, reorder_hook] = ismember([1:N_pruned_alleles-1]', reorder_hook);
    reorder_hook = [1; reorder_hook];

    [x, y, plot_order] = lineage_tree_layout(pruned_adj_list, pruned_allele_at_node, pruned_pop, ...
                                             reorder_hook, [1:max(reorder_hook)]');

    fig_width = 18.3;
    fig_height = 15;

    Nr = 2;
    Nc = 4;

    tight_margin = 0.2;
    top_margin = 0.4;
    bottom_margin = 0.8;

    pop_width = 1.0;
    pro_width = 4;
    tree_width = 4;
    bar_width = 3;
    allele_width = fig_width-pop_width-pro_width-tree_width-(Nc+1)*tight_margin;
    tissue_label_height = 1;
    matrix_height = pro_width;
    rel_height = fig_height-bottom_margin-top_margin-matrix_height-tissue_label_height;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
    
    bank1 = load(sprintf('%s/Banks/Protocol1/RNA/PosDox/Bank.mat', results_dir), 'bank'); bank1 = bank1.bank;
    bank2 = load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir), 'bank'); bank2 = bank2.bank;    
    [epval1, rate1] = bank1.compute_frequency_pvalue(raw_combined.summary);
    [epval2, rate2] = bank2.compute_frequency_pvalue(raw_combined.summary);
    epval = max(epval1, epval2);
    rate  = max(rate1,  rate2);
    
    sp1 = subplot(Nr, Nc, 1);
    scatter(log10(rate(2:end)*sum(raw_combined.summary.allele_freqs(2:end))), log10(raw_combined.summary.allele_freqs(2:end)), 6, log10(epval(2:end)), 'filled');    
    colormap(sp1, flipud(jet));
    cb6 = colorbar('FontSize', 5);    
    title(cb6, 'log_{10} pValue');
    caxis([-6 0]);
    axis tight;
    box on;
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('log_{10} Expected Frequency', 'FontSize', 6);
    ylabel('log_{10} Observed Frequency', 'FontSize', 6);
    title('Statistical Significance of Alleles', 'FontSize', 6', 'FontWeight', 'normal');
    
    sp2 = subplot(Nr, Nc, 2);
    temp = [sum(allele_breakdown_by_sample(2:end,:),1);
            sum(raw_combined.allele_breakdown_by_sample(2:end,:),1)]';
    barh(1:N_samples, [temp(:,1) diff(temp, [], 2)], 'stacked', 'BarWidth', 1.0);
    axis tight;    
    set(get(gca, 'XAxis'), 'Dir', 'reverse', 'FontSize', 5);
    set(get(gca, 'YAxis'), 'Dir', 'reverse', 'FontSize', 5);    
    xlabel('Edited Transcripts', 'FontSize', 6);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out', 'ytick', [], 'yticklabel', {});   
    lgd = legend({'Sig'; 'All'}, 'FontSize', 5);
    legend('boxoff');
    
    sp3 = subplot(Nr, Nc, 3);
    tissue_distance = compute_tissue_distance(stable_adj_list, allele_breakdown_by_sample(stable_allele_at_node(1:N_stable_alleles),:));
    imagesc(tissue_distance, 'XData', 1:N_samples, 'YData', 1:N_samples);
    axis tight; box on;
    colormap(sp3, 'pink');
    colorbar('FontSize', 5);
    set(gca, 'ydir', 'rev');
    set(gca, 'ytick', [1:N_samples], 'yticklabel', tissue_labels);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(gca, 'XAxisLocation', 'top', 'xtick', [], 'xticklabel', {});
    title('Tissue Distance', 'FontSize', 6, 'FontWeight', 'normal');
   
    sp4 = subplot(Nr,Nc,5);
    scatter(x, y, 3, 'black', 'filled'); 
    for i = 1:length(pruned_adj_list)
        if (pruned_adj_list(i) ~= 0)
            line([x(i) x(pruned_adj_list(i))], [y(i) y(pruned_adj_list(i))], 'color', [0.6 0.6 0.6]);
        end
    end
    set(gca, 'xdir', 'reverse');
    box off;
    axis off;
    axis tight;
       
    sp5 = subplot(Nr,Nc,6);
    RGB = get_sequence_coloring(flipud(aggregated_alleles(pruned_allele_at_node(plot_order(1:N_pruned_alleles)))), 'bp');
    imshow(padarray(RGB, [1 1], 0), 'YData', [0 N_pruned_alleles+1], 'XData', [0 CARLIN_def.getInstance.width.CARLIN+1]); 
    set(gca, 'ydir', 'normal');
    hold on;
    h = imshow(padarray(repmat([CARLIN_def.getInstance.alpha.CARLIN], [size(RGB,1), 1]), [1 1], 0), ...
               'YData', [0 N_pruned_alleles+1], 'XData', [0 CARLIN_def.getInstance.width.CARLIN+1]);  hold off;
    set(h, 'AlphaData', CARLIN_def.getInstance.alpha.overlay);
    axis tight;
    daspect([rel_height/(N_pruned_alleles+2), allele_width/(CARLIN_def.getInstance.width.CARLIN+2), 1]);
    
    sp6 = subplot(Nr,Nc,7);
    imagesc(padarray(tissue_frac(pruned_allele_at_node(plot_order(1:N_pruned_alleles)),:), [1, 0], 0), ...
            'XData', [1:N_samples], 'YData', [0:N_pruned_alleles+1]); colormap(sp6, pink);
    set(gca, 'ydir', 'normal');
    set(gca, 'ytick', [], 'yticklabel', {}, 'XAxisLocation', 'Top');
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(gca, 'xtick', [1:N_samples]);
    set(gca, 'xticklabels', tissue_labels);
    set(get(gca, 'XAxis'), 'TickLabelRotation', 90);
    axis tight; box off;
    cb = colorbar('Location', 'SouthOutside', 'FontSize', 5);    
    cb.Ticks = [0:4]/4;
    cb.AxisLocation = 'out';
        
    sp7 = subplot(Nr,Nc,8);
    pruned_pop_to_plot = log10(pruned_pop(plot_order(1:N_pruned_alleles)));
    barh(1:N_pruned_alleles, pruned_pop_to_plot, 'BarWidth', 1.0, 'EdgeColor', 'none', 'FaceColor', 'black');
    set(gca, 'ytick', [], 'yticklabel', {}, 'XAxisLocation', 'top');
    set(get(gca, 'XAxis'), 'FontSize', 5);
    axis tight; box off;
    xlim([0 ceil(max(pruned_pop_to_plot))]);
    set(gca, 'xtick', 0:ceil(max(pruned_pop_to_plot)));
        
    linkaxes([sp4, sp5, sp5, sp6, sp7], 'y');
    linkaxes([sp3, sp6], 'x');
    
    set(sp1, 'Units', 'centimeters', 'Position', [2*tight_margin+tree_width-0.8, bottom_margin+rel_height+tissue_label_height, pro_width, matrix_height]);
    set(sp2, 'Units', 'centimeters', 'Position', [2*tight_margin+tree_width-1+pro_width+1.5, bottom_margin+rel_height+tissue_label_height, bar_width, matrix_height]);    
    set(sp3, 'Units', 'centimeters', 'Position', [3*tight_margin+tree_width+allele_width, bottom_margin+rel_height+tissue_label_height, pro_width, matrix_height]);
    set(sp4, 'Units', 'centimeters', 'Position', [1*tight_margin bottom_margin, tree_width, rel_height]);
    set(sp5, 'Units', 'centimeters', 'Position', [2*tight_margin+tree_width, bottom_margin, allele_width, rel_height]);
    set(sp6, 'Units', 'centimeters', 'Position', [3*tight_margin+tree_width+allele_width, bottom_margin, pro_width, rel_height]);    
    set(sp7, 'Units', 'centimeters', 'Position', [4*tight_margin+tree_width+allele_width+pro_width, bottom_margin, pop_width, rel_height]);
        
    cb.Units = 'centimeters';
    cb.Position(4) = tight_margin;
    cb.Position(2) = 7*tight_margin/4;
    
    lgd.Units = 'centimeters';
    lgd.Position = [2*tight_margin+tree_width+pro_width+0.6, bottom_margin+rel_height+tissue_label_height+2, bar_width*0.6, matrix_height*0.1];    
        
    sample_label_width = 2.0;    
    sample_label = annotation('textbox', 'EdgeColor','none','string', 'Fraction In Tissue', ...
                              'Margin', 3, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 6);   
    sample_label.Units = 'centimeters';
    sample_label.Position = [3*tight_margin+tree_width+allele_width-sample_label_width, 0, sample_label_width, bottom_margin];
    
    transcript_label = annotation('textbox', 'EdgeColor','none','string', 'log_{10} Transcripts', ...
                                  'Margin', 0, 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'right', 'FontSize', 6);   
    transcript_label.Units = 'centimeters';
    transcript_label.Position = [fig_width-tight_margin-pop_width, 0 pop_width, bottom_margin];

    paper_print(sprintf('%s/ConsensusTree', outdir));
   
    plot_subtree_stability(depthX_subtree_instances{1}, depthX_subtree_clade_pop{1}, 1, N_sim);
    paper_print(sprintf('%s/Depth1Stability', outdir));
    
    plot_subtree_stability(depthX_subtree_instances{2}, depthX_subtree_clade_pop{2}, 2, N_sim);
    paper_print(sprintf('%s/Depth2Stability', outdir));
    
    plot_subtree_stability(depthX_subtree_instances{3}, depthX_subtree_clade_pop{3},3, N_sim);
    paper_print(sprintf('%s/Depth3Stability', outdir));
    
    close all;
    
end