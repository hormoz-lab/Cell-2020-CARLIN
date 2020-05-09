function stats = ks_comparison(pos_datavec, neg_datavec, N_ds, N_trials)

    [~, stats.pos.allele_pop] = resample_alleles(pos_datavec, N_ds, N_trials);        
    stats.pos.CDF = cellfun(@(x) cumsum(x)/sum(x), stats.pos.allele_pop, 'un', false);
    
    [~, stats.neg.allele_pop] = resample_alleles(neg_datavec, N_ds, N_trials);    
    stats.neg.CDF = cellfun(@(x) cumsum(x)/sum(x), stats.neg.allele_pop, 'un', false);
    
    stats.D_stat = max(abs(horzcat(stats.pos.CDF{:})-horzcat(stats.neg.CDF{:})), [], 1);
    stats.pval = percentile(exp(-N_ds*stats.D_stat.^2), 99);
    
end

