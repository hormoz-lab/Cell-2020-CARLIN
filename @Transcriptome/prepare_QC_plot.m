function prepare_QC_plot(counts, cell_mask, mt_genes, cutoff)

    keySet = {'Passed UMI/Passed MT'; 'Passed UMI/Failed MT'; 'Failed UMI/Passed MT';  'Failed UMI/Failed MT'};
    valSet = {'Blue'; 'Red'; 'Green'; 'Yellow'};    
    map_color = containers.Map(keySet, valSet);

    total_counts = full(sum(counts,2));
    mt_counts = full(sum(counts(:,mt_genes),2));

    figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot(2,3,1);
    [x, y] = split_series(log10(total_counts), cell_mask.umi, cell_mask.mt, 50);
    b = bar(x, y, 'stacked');
    [b.FaceColor] = valSet{:};
    line([log10(cutoff.umi), log10(cutoff.umi)], get(gca, 'ylim'), 'Color', map_color('Failed UMI/Passed MT'));
    title('UMI distribution'); ylabel('Number of Cells'); xlabel('log10(Number of UMIs)');    
    
    subplot(2,3,2);
    [x, y] = split_series(log10(full(sum(logical(counts),2))), cell_mask.umi, cell_mask.mt, 50);
    b = bar(x, y, 'stacked');
    [b.FaceColor] = valSet{:};
    title('Gene distribution'); ylabel('Number of Cells'); xlabel('log10(Number of Genes)');

    subplot(2,3,3);
    [x, y] = split_series(mt_counts./total_counts, cell_mask.umi, cell_mask.mt, 100);
    b = bar(x, y, 'stacked');
    [b.FaceColor] = valSet{:};
    line([cutoff.mt_frac, cutoff.mt_frac], get(gca, 'ylim'), 'Color', map_color('Passed UMI/Failed MT'));
    title('Mitochondrial UMI distribution'); ylabel('Number of Cells'); xlabel('Fraction of UMIs in MT genes');
    
    subplot(2,3,4);
    x = log10(total_counts);
    y = log10(full(sum(logical(counts),2)));    
    scatter(x( cell_mask.umi &  cell_mask.mt), y( cell_mask.umi &  cell_mask.mt), 6, map_color('Passed UMI/Passed MT'), 'filled'); hold on;
    scatter(x( cell_mask.umi & ~cell_mask.mt), y( cell_mask.umi & ~cell_mask.mt), 6, map_color('Passed UMI/Failed MT'), 'filled');
    scatter(x(~cell_mask.umi &  cell_mask.mt), y(~cell_mask.umi &  cell_mask.mt), 6, map_color('Failed UMI/Passed MT'), 'filled');
    scatter(x(~cell_mask.umi & ~cell_mask.mt), y(~cell_mask.umi & ~cell_mask.mt), 6, map_color('Failed UMI/Failed MT'), 'filled'); hold off;
    title('Gene vs UMI Distribution'); xlabel('log10(Number of UMIs)'); ylabel('log10(Number of Genes)');
    
    subplot(2,3,5);
    [x, y] = split_series(log10(full(sum(logical(counts),2)))./log10(total_counts), cell_mask.umi, cell_mask.mt, 50);
    b = bar(x, y, 'stacked');
    [b.FaceColor] = valSet{:};
    title('Gene per UMI Distribution'); xlabel('log10(Number of Genes)/log10(Number of UMIs)'); ylabel('Number of Cells');
    
    subplot(2,3,6);
    axis off;
    box off;
    text(0, 0.7, [sprintf('%-20s (%10s):\t%5d cells\n', keySet{1}, valSet{1}, sum(cell_mask.umi & cell_mask.mt)) ...
                  sprintf('%-20s (%10s):\t%5d cells\n', keySet{2}, valSet{2}, sum(cell_mask.umi & ~cell_mask.mt)) ...
                  sprintf('%-20s (%10s):\t%5d cells\n', keySet{3}, valSet{3}, sum(~cell_mask.umi & cell_mask.mt)) ...
                  sprintf('%-20s (%10s):\t%5d cells\n', keySet{4}, valSet{4}, sum(~cell_mask.umi & ~cell_mask.mt)) ...
                  sprintf('%-33s: %5d cells\n', 'Total', size(counts,1))], 'FontName', 'FixedWidth');
              
end           
              
function [edges, counts] = split_series(y, mask1, mask2, bins)

    assert(size(y,1) == size(mask1,1));
    assert(size(y,1) == size(mask2,1));
    
    y(y==-inf) = min(y(isfinite(y)));
    y(isnan(y)) = min(y(isfinite(y)));
    
    [~,edges] = histcounts(y, bins);
    counts = zeros(size(edges,2)-1,4);
    counts(:,1) = histcounts(y(mask1 & mask2), edges);
    counts(:,2) = histcounts(y(mask1 & ~mask2), edges);
    counts(:,3) = histcounts(y(~mask1 & mask2), edges);
    counts(:,4) = histcounts(y(~mask1 & ~mask2), edges);
    
    edges = movmean(edges, 2, 'Endpoints', 'discard')';
    
    assert(sum(counts(:)) == size(y,1));
end
