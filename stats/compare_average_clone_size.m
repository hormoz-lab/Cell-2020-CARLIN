function stats = compare_average_clone_size(pos_datavec, neg_datavec, N_trials)

    % We're computing average clone sizes through bootstrapping instead of
    % directly doing a t-test, to rule out the possibility that large clone
    % sizes are enabled purely through being able to sample more cells.
    % Since the 5-FU mice are characterized by a few large clones, bootstrapping
    % should flatten this distribution. Conversely, the controls are
    % relatively flat the sampling distribution should look more skewed.
    % Taken together this should yield a more conservative estimate.

    N_pos = length(pos_datavec);
    N_neg = length(neg_datavec);
   
    stats.pos.ds = resample_alleles(pos_datavec, N_pos, N_trials);
    stats.neg.ds = resample_alleles(neg_datavec, N_neg, N_trials);
  
    stats.pos.acs = N_pos./cellfun(@(x) length(unique(x)), stats.pos.ds');
    stats.neg.acs = N_neg./cellfun(@(x) length(unique(x)), stats.neg.ds');
    
    stats.pos.mean = length(pos_datavec)/length(unique(pos_datavec));
    stats.neg.mean = length(neg_datavec)/length(unique(neg_datavec));

    % To be conservative, if the mean due to bootstrapping is larger for
    % the positive treatment, shift the distribution to be the mean of the
    % original data. Likewise, if the mean for the bootstrapped version is
    % smaller for the negative treatment, shift it to the mean of the
    % original data.   
    
    if (mean(stats.pos.acs) > stats.pos.mean)
        stats.pos.acs = stats.pos.acs-mean(stats.pos.acs)+stats.pos.mean;
    end
    
    if (mean(stats.neg.acs) < stats.neg.mean)        
        stats.neg.acs = stats.neg.acs-mean(stats.neg.acs)+stats.neg.mean;
    end
    
    [stats.pval, stats.t, stats.df] = ttest_from_bootstrap(stats.pos.acs, stats.neg.acs, 'right', N_pos, N_neg);

end
