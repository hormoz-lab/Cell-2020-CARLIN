function make_invitro_phylogeny_subplots(results_dir)

    load(sprintf('%s/Trees/InVitroPhylogeny/Lineage.mat', results_dir));
    load(sprintf('%s/Trees/InVitroPhylogeny/Analysis.mat', results_dir));
    
    outdir = sprintf('%s/Figures/InVitroPhylogeny', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    [pruned_adj_list, pruned_allele_at_node, N_pruned_alleles] ...
        = generate_adj_list_from_subtrees(stable_depthX_subtrees, depthX_subtree_clade_pop, 'min', 100);
    
    % Just for making the plot, we'll show the ground truth library, as the
    % one in which an allele appears most if it appears in multiple libraries. Doesn't affect statistics.
    [~, ground_truth] = max(allele_breakdown_by_sample, [], 2);
    ground_truth = ground_truth+1;
    ground_truth(1) = 1;

    pruned_pop = sum(allele_breakdown_by_sample(pruned_allele_at_node(1:N_pruned_alleles),:),2);
    [x, y, plot_order] = lineage_tree_layout(pruned_adj_list, pruned_allele_at_node, pruned_pop, ...
                                             ground_truth(pruned_allele_at_node(1:N_pruned_alleles)), [1:max(ground_truth)]');
    
    [is_sanger, which_sanger] = ismember(pruned_allele_at_node, sanger_allele);

    fig_width  = 10.8;
    fig_height = 10.8;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    Nr = 1;
    Nc = 3;

    tight_margin = 0.1;
    bottom_margin = 0.4;

    pop_width = 1.0;
    tree_width = 3;
    allele_width = fig_width-pop_width-tree_width-(Nc+1)*tight_margin;
    rel_height = fig_height-bottom_margin-tight_margin;

    % Match colors from schematic
    color_scheme = [128 128 128;
                    233 185 210;
                    240 185 140;
                    163 181 141;
                     83 177 247;
                    169 145 145;
                    212  94  94;
                    152 174 246;
                     91 165 100]/255;

    sp1 = subplot(Nr,Nc,1);
    
    hold on;
    colors = color_scheme(ground_truth(pruned_allele_at_node),:);
    for i = 1:length(pruned_allele_at_node)
        if (pruned_adj_list(i) ~= 0)
            line([x(i) x(pruned_adj_list(i))], [y(i) y(pruned_adj_list(i))], 'color', colors(i,:));
        end
    end    
    scatter(x, y, 3, colors, 'filled'); 
    scatter(x(is_sanger), y(is_sanger), 6, color_scheme(which_sanger(is_sanger),:), 'filled'); 

    set(gca, 'xdir', 'reverse');
    box off;
    axis off;
    axis tight;
    hold off;
    
    sp2 = subplot(Nr,Nc,2);
    RGB = get_sequence_coloring(flipud(aggregated_alleles(pruned_allele_at_node(plot_order(1:N_pruned_alleles)))), 'bp');
    imshow(padarray(RGB, [1 1], 0), 'YData', [0 N_pruned_alleles+1], 'XData', [0 CARLIN_def.getInstance.width.CARLIN+1]);
    set(gca, 'ydir', 'normal');
    hold on;
    h = imshow(padarray(repmat([CARLIN_def.getInstance.alpha.CARLIN], [size(RGB,1), 1]), [1 1], 0), ...
               'YData', [0 N_pruned_alleles+1], 'XData', [0 CARLIN_def.getInstance.width.CARLIN+1]);  hold off;
    set(h, 'AlphaData', CARLIN_def.getInstance.alpha.overlay);
    axis tight;
    daspect([fig_height*rel_height/(N_pruned_alleles+2), fig_width*allele_width/(CARLIN_def.getInstance.width.CARLIN+2), 1]);

    
    sp3 = subplot(Nr,Nc,3);
    barh(log10(pruned_pop(plot_order(1:N_pruned_alleles))), 'BarWidth', 1.0, 'EdgeColor', 'none', 'FaceColor', 'black');
    hold on;    
    for i = find(is_sanger(1:N_pruned_alleles))'    
        barh(find(plot_order==i), log10(pruned_pop(i)), 'BarWidth', 1.0, 'EdgeColor', 'none', 'FaceColor', color_scheme(which_sanger(i),:));
    end
    hold off;
    set(gca, 'ytick', [], 'yticklabel', {}, 'XAxisLocation', 'bottom');
    set(get(gca, 'XAxis'), 'FontSize', 5);
    axis tight; box off;
    maxtick = ceil(log10(max(pruned_pop)));
    xlim([0 maxtick]);    
    set(gca, 'xtick', [0:maxtick], 'xticklabel', cellstr(num2str([0:maxtick]')));

    linkaxes([sp1, sp2, sp3], 'y');
    
    set(sp1, 'Units', 'centimeters', 'Position', [tight_margin, bottom_margin, tree_width, rel_height]);
    set(sp2, 'Units', 'centimeters', 'Position', [2*tight_margin+tree_width, bottom_margin, allele_width, rel_height]);
    set(sp3, 'Units', 'centimeters', 'Position', [3*tight_margin+tree_width+allele_width, bottom_margin, pop_width, rel_height]);
    
    label_width = 1.1;
    
    an = annotation('textbox', 'EdgeColor','none','string', 'log_{10} Cells', ...
                    'Margin', 0, 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'left', 'FontSize', 6);   
    an.Units = 'centimeters';
    an.Position = [fig_width-tight_margin-pop_width-tight_margin-label_width, 0, label_width, bottom_margin];

    paper_print(sprintf('%s/ConsensusTree', outdir));
    
    fig_width = 3;
    fig_height = 4;
    left_margin = 1.3;
    bottom_margin = 0.5;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
    
    sp = subplot(1,1,1);
    imagesc(master_event_table, 'XData', 1:N_samples, 'YData', 1:size(master_sanger_events,1)); 
    colormap([1 1 1; color_scheme(2:end,:)]);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(gca, 'XAxisLocation', 'top');
    set(gca,'TickLabelInterpreter','none');
    xticks(1:N_samples); 
    xlabel('Sanger Allele', 'FontSize', 6);
    yticks(1:size(master_sanger_events,1));
    yticklabels(master_sanger_events);   
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin 0 fig_width-left_margin fig_height-bottom_margin]);
    paper_print(sprintf('%s/SangerMutationTable', outdir));    
   
    plot_subtree_stability(depthX_subtree_instances{1}, depthX_subtree_clade_pop{1}, 1, N_sim);
    paper_print(sprintf('%s/Depth1Stability', outdir));
    
    plot_subtree_stability(depthX_subtree_instances{2}, depthX_subtree_clade_pop{2}, 2, N_sim);
    paper_print(sprintf('%s/Depth2Stability', outdir));
    
    plot_subtree_stability(depthX_subtree_instances{3}, depthX_subtree_clade_pop{3},3, N_sim);
    paper_print(sprintf('%s/Depth3Stability', outdir));
    
    close all;
    
end
