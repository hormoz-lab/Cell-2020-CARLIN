function [KS_stat, KS_dist] = bootstrap_KS_statistic(counts, N_trials)
    
    assert(size(counts,2) == 2);
    N_bins = size(counts,1);
    N_obs = sum(counts,1);
    
    KS_stat = max(abs(diff(cumsum(counts,1)./N_obs, [], 2)));
    
    datavec1 = repelem([1:N_bins]', counts(:,1));
    datavec2 = repelem([1:N_bins]', counts(:,2));
    
    rs_counts{1} = arrayfun(@(i) accumarray(datasample(datavec1, N_obs(1), 'Replace', true), 1, [N_bins, 1]), [1:N_trials]', 'un', false);
    rs_counts{2} = arrayfun(@(i) accumarray(datasample(datavec2, N_obs(2), 'Replace', true), 1, [N_bins, 1]), [1:N_trials]', 'un', false);
    
    KS_dist = sort(cellfun(@(x,y) max(abs(diff(cumsum([x,y],1)./N_obs, [], 2))), rs_counts{1}, rs_counts{2}));

end