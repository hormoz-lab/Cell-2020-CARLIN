function [umi_mask, umi_cutoff] = filter_out_cells_low_umi(A, umi_cutoff)

    fprintf('Filtering low UMI cells: ');

    total_counts = sum(A,2);
    
    if (nargin < 2 || umi_cutoff < 0)
        if (min(total_counts) >= 1000)
            umi_cutoff = min(total_counts);
        else
            [h, b] = histcounts(log10(total_counts), 50);
            b = movmean(b,2, 'Endpoints', 'discard')';
            [~, loc, width, prom] = findpeaks(-h');
            sorted = sortrows([prom, width, b(loc)], [1, 2, 3], {'Descend'; 'Ascend'; 'Ascend'});
            umi_cutoff = ceil(10^sorted(find(sorted(:,3)>=3.0 & sorted(:,3) < 3.5, 1, 'first'),3));
        end
        umi_cutoff = full(umi_cutoff);
    end
    
    umi_mask = full(total_counts >= umi_cutoff);
    fprintf('%d/%d survive at UMI cutoff >= %d\n', sum(umi_mask), length(umi_mask), umi_cutoff);
    
end
    
    
    
    