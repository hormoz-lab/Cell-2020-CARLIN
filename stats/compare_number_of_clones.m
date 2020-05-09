function stats = compare_number_of_clones(pos, neg)

    assert(isa(pos.summary, 'ExperimentSummary'));
    assert(isa(neg.summary, 'ExperimentSummary'));
    
    stats.pos.datavec = repelem([1:length(pos.sig_alleles.all)]', pos.summary.allele_freqs(pos.sig_alleles.all));
    stats.neg.datavec = repelem([1:length(neg.sig_alleles.all)]', neg.summary.allele_freqs(neg.sig_alleles.all));
    
    stats.N_trials = 10000;
    
    stats.N_ds = min(length(stats.pos.datavec), length(stats.neg.datavec));
    
    [~, ~, stats.pos.csf] = resample_alleles(stats.pos.datavec, stats.N_ds, stats.N_trials);
    [~, ~, stats.neg.csf] = resample_alleles(stats.neg.datavec, stats.N_ds, stats.N_trials);
    
    stats.pos.nc = cellfun(@sum, stats.pos.csf);
    stats.neg.nc = cellfun(@sum, stats.neg.csf);
    
    % We downsample to control the effect of seeing fewer clones when there
    % are fewer cells. However, bootstrapping only samples on average 1/e 
    % of the samples on each iteration, which can artificially lower the number
    % of detected clones. A correction for this is to re-shift the distribution
    % to the number of detected in the original sample. If the smaller sample 
    % is the positive treatment, this will inflate the number of clones, which is conservative. 
    % If it's the negative treatment, we won't apply this correction, so that the
    % number of clones is under-estimated, which is again conservative.
    
    if (stats.N_ds == length(stats.pos.datavec))
        stats.pos.nc = stats.pos.nc-mean(stats.pos.nc)+length(unique(stats.pos.datavec));  
    end
    
    [stats.pval, stats.t, stats.df] = ttest_from_bootstrap(stats.pos.nc, stats.neg.nc, 'left', stats.N_ds);
   
end