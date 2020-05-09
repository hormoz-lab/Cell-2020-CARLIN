function make_diversity_subplots(results_dir)
    
    outdir = sprintf('%s/Figures/Diversity', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    load(sprintf('%s/Banks/Protocol2/RNA/NegDox/Bank.mat', results_dir));

    plot_highlighted_alleles(bank.summary, 50);
    paper_print(sprintf('%s/NegDoxTop50Sequences', outdir));
    
    plot_allele_frequency_CDF(bank.summary, '-Dox');
    paper_print(sprintf('%s/NegDoxCDF', outdir));
    
    load(sprintf('%s/Banks/Protocol1/RNA/NegCas9/Bank.mat', results_dir));
    
    plot_highlighted_alleles(bank.summary, length(bank.summary.alleles)-1);
    paper_print(sprintf('%s/NegCas9Top50Sequences', outdir));
    
    plot_allele_frequency_CDF(bank.summary, '-Cas9');
    paper_print(sprintf('%s/NegCas9CDF', outdir));
        
    load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir));
    
    plot_highlighted_alleles(bank.summary, 50);
    paper_print(sprintf('%s/PosDoxTop50Sequences', outdir));
    
    plot_allele_frequency_CDF(bank.summary, '+Dox');
    paper_print(sprintf('%s/PosDoxCDF', outdir));
    
    plot_indel_freq_vs_length(bank.summary, bank.allele_breakdown_by_sample);
    paper_print(sprintf('%s/BankInDelFreqVsLength', outdir));
    
    plot_allele_membership_overlap(bank);
    paper_print(sprintf('%s/BankAlleleMembership', outdir));
    
    plot_theoretical_allele_discovery_curve(bank);    
    paper_print(sprintf('%s/BankTheoreticalExtrapolation', outdir));
    
    rng(10,'twister');    
    plot_empirical_allele_discovery_curve(bank);
    paper_print(sprintf('%s/BankEmpiricalExtrapolation', outdir));

    M = [500:500:50000]';
    sig_level = 0.05;
    N = length(M);
    alleles_total_expected = zeros(N,1);
    rate_at_sig_level = zeros(N,1);
    alleles_sig_expected = zeros(N,1);

    for i = 1:length(M)
        alleles_total_expected(i) = bank.interpolate_alleles(M(i));
        rate_at_sig_level(i) = cell2mat(Bank.max_rate_for_sig_level(sig_level, M(i)));    
        alleles_sig_expected(i) = alleles_total_expected(i)*integral(@(x) (1-exp(-x*M(i))).*bank.model.intensive.rate_pdf(x), 0, rate_at_sig_level(i))/ ...
                                                            integral(@(x) (1-exp(-x*M(i))).*bank.model.intensive.rate_pdf(x), 0, bank.model.intensive.cutoff);
    end
    syms M_clamp;
    clamp_val = floor(double(vpasolve(Bank.clonal_pvalue_eqn(bank.model.intensive.mean_unobs_rate, M_clamp)==sig_level)));
    plot_useful_alleles_calibration(M, alleles_sig_expected, clamp_val);
    paper_print(sprintf('%s/InductionCalibration', outdir));
    
    
    load(sprintf('%s/Simulation/DiversityReplicate/Results.mat', results_dir), ...
        'N_sampled', 'N_unique1', 'N_unique2', 'N_shared', 'N_union', 'pval_union');    
    plot_diversity_replicates(N_sampled, N_unique1, N_unique2, N_shared, N_union, pval_union);    
    paper_print(sprintf('%s/ReplicatePValue', outdir));
    
    close all;
end