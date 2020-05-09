function stats = compare_clone_size_distribution(pos, neg)

    assert(isa(pos.summary, 'ExperimentSummary'));
    allele_freqs = pos.summary.allele_freqs;
    pos = compute_hsc_clone_sizes(pos);
    pos.datavec = repelem([1:length(pos.test_alleles)]', allele_freqs(pos.test_alleles));
    
    stats.N_trials = 10000;
    
    % If we're provided two samples, compare average clone size between the
    % two using a t-test
    if (nargin == 2)
        assert(isa(neg.summary, 'ExperimentSummary'));
        allele_freqs = neg.summary.allele_freqs;
        neg = compute_hsc_clone_sizes(neg);
        neg.datavec = repelem([1:length(neg.test_alleles)]', allele_freqs(neg.test_alleles));
        stats = merge_structs(stats, compare_average_clone_size(pos.datavec, neg.datavec, stats.N_trials));       
        
    % Or else compare the sample to a uniform using a KS-test
    else
        neg.datavec = [1:length(pos.test_alleles)]';
        stats.N_ds = pos.N_transcripts;
        stats = merge_structs(stats, ks_comparison(pos.datavec, neg.datavec, stats.N_ds, stats.N_trials));        
    end   
    
     stats.pos = merge_structs(pos, stats.pos);
     stats.neg = merge_structs(neg, stats.neg);
end