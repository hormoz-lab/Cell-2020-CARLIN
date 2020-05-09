function plot_number_clones_downsampled(samples, colors, xlabel_str, title_str)

    assert(length(samples)==2);
    assert(length(xlabel_str)==2);

    which_clones = cellfun(@(x) x.sig_alleles.all, samples, 'un', false);
    clone_size   = cellfun(@(x,y) x.summary.allele_freqs(y), samples, which_clones, 'un', false);
    total_size  = cellfun(@(x) sum(x), clone_size);
    n = min(total_size);
    ds = find(total_size~=n);
        
    rng(39594);
    which_clones{ds} = datasample(repelem(which_clones{ds}, clone_size{ds}), n, 'replace', false);
    [which_clones{ds}, ~, clone_size{ds}] = find(accumarray(which_clones{ds}, 1));
    assert(sum(clone_size{ds})==n);

    pad_to_size = max(length(clone_size{1}), length(clone_size{2})); 
    which_clones = cellfun(@(x) padarray(x, pad_to_size-length(x), NaN, 'post'), which_clones, 'un', false);
    clone_size   = cellfun(@(x) padarray(x, pad_to_size-length(x), NaN, 'post'), clone_size, 'un', false);
    clone_type = cell(2,1);
    
    for j = 1:2
        clone_type{j} = NaN(size(which_clones{j}));
        for i = 1:length(which_clones{j})
            if (ismember(which_clones{j}(i), samples{j}.hsc_derived_allele))
                clone_type{j}(i) = 1;            
            elseif (ismember(which_clones{j}(i), samples{j}.hsc_only_allele))
                clone_type{j}(i) = 2;
            elseif (ismember(which_clones{j}(i), samples{j}.derived_only_allele))
                clone_type{j}(i) = 3;
            end
        end
        [~, reorder] = sortrows([clone_size{j} clone_type{j}], [1 2], {'descend'; 'ascend'});
        clone_size{j} = clone_size{j}(reorder);
        clone_type{j} = clone_type{j}(reorder);
        which_clones{j} = which_clones{j}(reorder);
    end
    
    fig_width  = 1.5;
    fig_height = 4.2;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
  
    hold on;
    for j = 1:2
        b{j} = bar([j; NaN], [clone_size{j} NaN(size(clone_size{j}))]', 'stacked', 'BarWidth', 1, 'LineWidth', 0.1);
        for i = 1:length(b{j})
            if (~isnan(clone_type{j}(i)))            
                set(b{j}(i),'FaceColor', colors(clone_type{j}(i),:));
            end
        end
    end
    hold off;
    axis tight;
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(gca, 'Xtick', [1:2], 'xticklabel', xlabel_str, 'ytick', []);
    xtickangle(45);
    
    ylabel(sprintf('Fraction of Cells (n=%d)', n), 'FontSize', 6);
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    set(gca, 'LineWidth', 0.1);
    set(gca, 'TickDir', 'out');
    
    left_margin = 0.4;
    right_margin = 0.0;
    top_margin = 0.3;
    bottom_margin = 1.0;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end