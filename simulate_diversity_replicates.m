function simulate_diversity_replicates(analysis_dir, results_dir)

    load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir));

    SB213 = load(sprintf('%s/Diversity/SB213/Combined/Summary.mat', analysis_dir), 'summary');
    SB214 = load(sprintf('%s/Diversity/SB214/Combined/Summary.mat', analysis_dir), 'summary');
    [combined, ~, allele_breakdown_by_sample] = ExperimentSummary.FromMerge([SB213.summary; SB214.summary]);
    
    [~, rates] = bank.compute_clonal_pvalue(combined);
    datavec1 = repelem([2:length(combined.alleles)]', allele_breakdown_by_sample(2:end,1));
    datavec2 = repelem([2:length(combined.alleles)]', allele_breakdown_by_sample(2:end,2));
    N_sampled = [100:100:1000 2000:1000:10000 20000:10000:50000]';
    
    rng(3023940);
    
    N_trials = 1000;
    ds1 = arrayfun(@(x) resample_alleles(datavec1, x, N_trials), N_sampled/2, 'un', false);
    ds2 = arrayfun(@(x) resample_alleles(datavec2, x, N_trials), N_sampled/2, 'un', false);
    ds1 = vertcat(ds1{:});
    ds2 = vertcat(ds2{:});
    ds1 = cellfun(@(x) accumarray(x,1), ds1, 'UniformOutput', false);
    ds2 = cellfun(@(x) accumarray(x,1), ds2, 'UniformOutput', false);
    
    N_unique1 = cellfun(@(x,y) setdiff(find(x), find(y)), ds1, ds2, 'un', false);
    N_unique2 = cellfun(@(x,y) setdiff(find(y), find(x)), ds1, ds2, 'un', false);
    N_shared  = cellfun(@(x,y) intersect(find(x), find(y)), ds1, ds2, 'un', false);
    N_union  = cellfun(@(x,y) union(find(x), find(y)), ds1, ds2, 'un', false);
    pval_union = cellfun(@(r,M) bank.clonal_pvalue_eqn(rates(r), M), N_union, num2cell(repmat(N_sampled, [1, N_trials])), 'UniformOutput', false);
    
    clearvars -except N_sampled N_unique1 N_unique2 N_shared N_union pval_union results_dir
    
    outdir = sprintf('%s/Simulation/DiversityReplicate', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    save(sprintf('%s/Results.mat', outdir), 'N_sampled', 'N_unique1', 'N_unique2', 'N_shared', 'N_union', 'pval_union');
    
end
