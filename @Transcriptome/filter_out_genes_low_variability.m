function [masks, vscores] = filter_out_genes_low_variability(A, min_mean, min_count, min_cell, v_pctl)

    fprintf('Filtering low variability genes\n');
    
    N = size(A,2);

    % 1. Only compute vscores on genes that have exceed a given mean expression
    % level
    mean_gene = full(mean(A,1))';
    masks.exceeds_min_mean = mean_gene > min_mean;
    
    % 2. FF on genes that are at least minimally expressed
    
    mean_gene = mean_gene(masks.exceeds_min_mean);        
    var_gene = full(var(A(:,masks.exceeds_min_mean),0,1))';    
    FF_gene = var_gene./mean_gene;
    
    % 3. Compute variability as a function of mean expression
    
    sorted = sortrows([log(mean_gene), log(FF_gene./mean_gene)], [1 2], {'Ascend', 'Ascend'});
    
    bin_edges = linspace(min(sorted(:,1)), max(sorted(:,1)), 51)';
    bin_lo = bin_edges(1:end-1);
    bin_hi = bin_edges(2:end);
    
    x = (bin_lo+bin_hi)/2.0;
    y = arrayfun(@(i) percentile(sorted(sorted(:,1) >= bin_lo(i) & sorted(:,1) < bin_hi(i),2), 0.10), [1:size(x,1)]');
    nan_mask = isnan(y);
    y(nan_mask) = y(arrayfun(@(i) find(~nan_mask(1:i), 1, 'last'), find(nan_mask)));
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot(2,2,1);
    scatter(sorted(:,1), sorted(:,2), '.'); hold on;
    plot(x,y); hold off;
    title('(FF/Mean) vs (Mean) for genes');
    ylabel('log(FF/Mean)'); xlabel('log(Mean)');
    legend({'All genes'; '10th %-ile for bin'});
    
    % 4. FF distribution for genes
    
    [h,bin_edges] = histcounts(log(FF_gene), 200);    
    bin_values = movmean(bin_edges,2,'Endpoints', 'discard');
    
    subplot(2,2,2);
    bar(bin_values, h);
    title('FF distribution for genes'); ylabel('Frequency'); xlabel('log(FF)');
    
    % 5. Fit percentile curve for FF vs mean
    
    c = max(exp(bin_values(h==max(h))),1);
    
    gLog = @(x,c,b) log(c*exp(-x)+b);
    err = @(b2) sum(abs(gLog(x,c,b2)-y));
    options = optimset('Display','iter');
    b = fminsearch(err, 0.1, options);
    a = c / (1+b)-1;
    
    vscores = FF_gene ./ (c + b * mean_gene);
    CV_eff = sqrt(c - 1);
    CV_input = sqrt(b);
    
    masks.nonzero_vscore = false(N,1);
    masks.nonzero_vscore(masks.exceeds_min_mean) = (vscores > 0);

    mean_gene = mean_gene(vscores>0);
    FF_gene = FF_gene(vscores>0);
    vscores = vscores(vscores>0);        
    min_vscore = percentile(vscores, v_pctl);
    
    masks.vscore_exceeds_min = false(N,1);
    masks.vscore_exceeds_min(masks.nonzero_vscore) = vscores' >= min_vscore;
    
    temp = vscores;
    vscores = zeros(N,1);
    vscores(masks.nonzero_vscore) = temp;
    clear temp;
    
    masks.exceeds_min_count = (sum(A >= min_count,1) >= min_cell)';
    masks.hv_genes = masks.exceeds_min_mean & masks.nonzero_vscore & masks.vscore_exceeds_min & masks.exceeds_min_count;
        
    x_min = min(mean_gene)/2;
    x_max = max(mean_gene)*2;
    xTh = x_min * exp(log(x_max/x_min)*linspace(0,1,100));
    yTh = c + b * xTh;
    
    subselect = masks.vscore_exceeds_min(masks.nonzero_vscore) & masks.exceeds_min_count(masks.nonzero_vscore);
    assert(sum(subselect) == sum(masks.hv_genes));
    
    subplot(2,2,3);
    scatter(log10(mean_gene), log10(FF_gene), '.'); hold on;
    scatter(log10(mean_gene(subselect)), log10(FF_gene(subselect)), '.');    
    plot(log10(xTh), log10(yTh));
    hold off;
    legend({'Non-zero Vscore'; sprintf('Vscore >= %d', v_pctl); 'Fitted'});    
    title('FF vs Mean for genes with non-zero Vscore');
    xlabel('log10(Mean)'); ylabel('log10(FF)');
    
    masks.exceeds_min_mean   = uint16(find(masks.exceeds_min_mean));
    masks.exceeds_min_count  = uint16(find(masks.exceeds_min_count));
    masks.nonzero_vscore     = uint16(find(masks.nonzero_vscore));
    masks.vscore_exceeds_min = uint16(find(masks.vscore_exceeds_min));
    masks.hv_genes           = uint16(find(masks.hv_genes));
    
    subplot(2,2,4);
    axis off;
    box off;
    text(0, 0.7, [sprintf('%-30s:\t%7d genes\n', 'All', N) ...
                  sprintf('%-30s:\t%7d genes\n', sprintf('Mean UMI > %d', min_mean), length(masks.exceeds_min_mean)) ...
                  sprintf('%-30s:\t%7d genes\n', sprintf('UMI>=%d in >= %d cells', min_count, min_cell), length(masks.exceeds_min_count)) ...
                  sprintf('%-30s:\t%7d genes\n', 'Vscore > 0', length(masks.nonzero_vscore)) ...
                  sprintf('%-30s:\t%7d genes\n', sprintf('Vscore %%-ile >= %4.2g', v_pctl), length(masks.vscore_exceeds_min)) ...
                  sprintf('%-30s:\t%7d genes\n', 'Highly variable', length(masks.hv_genes))], 'FontName', 'FixedWidth');

end

