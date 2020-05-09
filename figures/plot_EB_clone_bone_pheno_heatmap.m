function plot_EB_clone_bone_pheno_heatmap(clone_by_pheno_by_bone_breakdown, alleles_to_plot, bb_alleles, bb_pval, ...
                                          hsc_derived_allele, hsc_only_allele, derived_only_allele, colors, ...
                                          coarse_grain_pheno, title_str)

    fig_width = 10.6;
    fig_height = 5;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    [N_alleles, N_tissues, N_bones] = size(clone_by_pheno_by_bone_breakdown);
    assert(length(alleles_to_plot) == N_alleles);
    assert(length(title_str) == N_bones);
    assert(length(coarse_grain_pheno) == N_tissues);
       
    top_margin = 0.3;
    bottom_margin = 0.6;
    right_margin = 0.4;
    left_margin = 0.3;
    
    pop_width = 2.0;
    label_width = 1.0;
    cb_width = 0.2;
    tight_margin = 0.1;
    
    heatmap_width = (fig_width-left_margin-right_margin-pop_width-label_width-cb_width-tight_margin)/N_bones;
    rel_height = fig_height-top_margin-bottom_margin;
    
    sp = cell(N_bones+1,1);
    
    [is, where] = ismember(bb_alleles, alleles_to_plot);
    assert(all(is));
    pval_label = ones(size(is));
    pval_label(bb_pval < 1e-3) = 2;
    pval_label(bb_pval < 1e-6) = 3;
    pval_label = arrayfun(@(i) repelem('*', i), pval_label, 'un', false);
    
    sp{1} = subplot(1,N_bones+1,1);
    b = barh(1:N_alleles, fliplr(squeeze(sum(clone_by_pheno_by_bone_breakdown,2))), 'stacked');
    ax = gca;
    set(ax, 'ydir', 'reverse', 'xdir', 'reverse', 'yaxislocation', 'right');
    axis tight;    
    set(ax, 'ytick', where, 'yticklabel', pval_label);
    hold on;
    box off;
    yyaxis right;
    set(ax, 'ydir', 'reverse', 'xdir', 'reverse');
    set(ax, 'ytick', [1:N_alleles]', 'yticklabel', arrayfun(@(x) sprintf('Clone %d', x), alleles_to_plot, 'un', false));
    xl = xlim;
    xlim([0, xl(2)])
    ax = gca;
    
    for j = 1:length(ax.YTickLabel)
        if (ismember(alleles_to_plot(j), hsc_derived_allele))
            color_prefix = sprintf('\\color[rgb]{%f,%f,%f}', colors(1,1), colors(1,2), colors(1,3));
        elseif (ismember(alleles_to_plot(j), hsc_only_allele))
            color_prefix = sprintf('\\color[rgb]{%f,%f,%f}', colors(2,1), colors(2,2), colors(2,3));            
        elseif (ismember(alleles_to_plot(j), derived_only_allele))
            color_prefix = sprintf('\\color[rgb]{%f,%f,%f}', colors(3,1), colors(3,2), colors(3,3));        
        end
        ax.YTickLabel{j} = [color_prefix ax.YTickLabel{j}];
    end

   
    set(get(ax, 'XAxis'), 'FontSize', 5);
    set(get(ax, 'YAxis'), 'FontSize', 5, 'Color', 'black');
    box off;
    linkprop([ax.YAxis(1) ax.YAxis(2)], 'Limits');
    
    title('Clone Size', 'FontSize', 6, 'FontWeight', 'normal');
    
    cmap = ones(64,3);
    cmap(:,2) = linspace(0, 1, 64);
    cmap(:,3) = linspace(0, 1, 64);
    colormap(flipud(cmap));
    max_count = max(clone_by_pheno_by_bone_breakdown(:));
    
    for i = 1:N_bones
        sp{i+1} = subplot(1,N_bones+1,i+1);        
        imagesc(squeeze(clone_by_pheno_by_bone_breakdown(:,:,i)), 'XData', 1:N_tissues, 'YData', 1:N_alleles);        
        caxis([0 max_count]);
        set(gca, 'YTick', [], 'Xtick', 1:N_tissues, 'xticklabel', coarse_grain_pheno, 'TickDir', 'out');
        set(get(gca, 'XAxis'), 'FontSize', 5);
        xtickangle(90);
        title(title_str{i}, 'FontSize', 6, 'FontWeight', 'normal', 'color', b(N_bones-i+1).FaceColor);        
    end
    
    cb = colorbar('Location', 'EastOutside', 'FontSize', 5);
    cb.AxisLocation = 'out';
    
    linkaxes(horzcat(sp{:}), 'y');
    
    set(sp{1}, 'Units', 'centimeters', 'Position', [left_margin bottom_margin pop_width rel_height]);
    
    for i = 1:N_bones
        set(sp{i+1}, 'Units', 'centimeters', 'Position', [left_margin+pop_width+label_width+(i-1)*heatmap_width bottom_margin heatmap_width rel_height]);        
    end
    
    cb.Units = 'centimeters';
    cb.Position(1) = left_margin+pop_width+label_width+N_bones*heatmap_width+tight_margin;
    cb.Position(2) = bottom_margin;
    cb.Position(3) = cb_width;
    cb.Position(4) = rel_height;
    title(cb, 'Cells', 'FontSize', 6, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');

end   