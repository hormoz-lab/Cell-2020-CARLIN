function make_pulsed_induction_subplots(analysis_dir, results_dir)

    outdir = sprintf('%s/Figures/PulsedInduction', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
 
    N = 4;
    
    dat = cell(8,1);
    dat{1} = load(sprintf('%s/PulsedInduction/G0/Summary.mat', analysis_dir), 'summary');
    dat{2} = load(sprintf('%s/PulsedInduction/G1/Summary.mat', analysis_dir), 'summary');
    dat{3} = load(sprintf('%s/PulsedInduction/G2A/Summary.mat', analysis_dir), 'summary');
    dat{4} = load(sprintf('%s/PulsedInduction/G2B/Summary.mat', analysis_dir), 'summary');
    dat{5} = load(sprintf('%s/PulsedInduction/G3A1/Summary.mat', analysis_dir), 'summary');
    dat{6} = load(sprintf('%s/PulsedInduction/G3A2/Summary.mat', analysis_dir), 'summary');
    dat{7} = load(sprintf('%s/PulsedInduction/G3B1/Summary.mat', analysis_dir), 'summary');
    dat{8} = load(sprintf('%s/PulsedInduction/G3B2/Summary.mat', analysis_dir), 'summary');
    
    dat = cellfun(@(x) x.summary, dat, 'un', false);
        
    chord = cell(N,1);
    
    chord{1} = ExperimentSummary.FromMerge(dat{1});
    chord{2} = ExperimentSummary.FromMerge(dat{2});
    chord{3} = ExperimentSummary.FromMerge([dat{3}; dat{4}]);
    chord{4} = ExperimentSummary.FromMerge([dat{5}; dat{6}; dat{7}; dat{8}]);
    clear dat;
  
    plot_quad_stargate(chord, {'R0'; 'R1'; 'R2'; 'R3'});
    paper_print(sprintf('%s/Stargate', outdir));
    
    muts = cellfun(@(y) cellfun(@(x) Mutation.identify_Cas9_events(x), y.alleles, 'un', false), chord(2:4), 'un', false);
    
    N_muts = cellfun(@(x) cellfun(@length, x), muts, 'un', false);
    max_mut = max(vertcat(N_muts{:}));
    N_mut_alleles = cellfun(@(x) accumarray(x+1,1, [max_mut+1 1]), N_muts, 'un', false);
    N_mut_alleles = horzcat(N_mut_alleles{:});
    N_mut_alleles = N_mut_alleles(2:end,:);
    N_mut_alleles = N_mut_alleles./sum(N_mut_alleles,1);
    plot_mutations_after_induction_round(N_mut_alleles, {'R1'; 'R2'; 'R3'});
    paper_print(sprintf('%s/MutationAccrualWithRound', outdir));
    
    ee = zeros(N,1);
    apc = struct('mu', zeros(N,1), 'sig', zeros(N,1));
    di = struct('mu', zeros(N,1), 'sig', zeros(N,1));
    CP = struct('mu', zeros(N,1), 'sig', zeros(N,1), 'dist', {cell(N,1)});
    num_bp_del = struct('mu', zeros(N,1), 'sig', zeros(N,1));
	num_bp_ins = struct('mu', zeros(N,1), 'sig', zeros(N,1));
    edited_only = true;
  
    N_bootstrap_trials = 1000;
    rng(23845498);

    for i = 1:N        
        summary = chord{i};
        ee(i) = sum(summary.allele_freqs(2:end))/sum(summary.allele_freqs)*100;
        [apc.mu(i), apc.sig(i)]               = bootstrap_alleles_per_cell(summary, N_bootstrap_trials);     
        [di.mu(i), di.sig(i)]                 = bootstrap_diversity_index(summary, N_bootstrap_trials);
        [CP.mu(i), CP.sig(i), CP.dist{i}] = CARLIN_potential(summary);
        [num_bp_del.mu(i), num_bp_ins.mu(i), num_bp_del.sig(i), num_bp_ins.sig(i)] = ...
            Mutation.num_bps_indel_stats(summary, edited_only);
    end
        
    dat.series_names = {'Induction Round'};
    dat.xvals = [{[0:3]}];
    dat.markers = {'o'};
    
    plt.xlim = [-0.5 3.5];    
    plt.xticklabels = {'R0'; 'R1'; 'R2'; 'R3'};
    plt.xlabel = 'Induction Round';
    
    plt.title = 'BPs Deleted In Edited Cells';
    plt.ylabel = 'Average BPs Deleted';    
    dat.yvals = {num_bp_del.mu};
    plt.ylim = [0 250];

    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/BPsDeletedOverTime', outdir));
    
    plt.title = 'Allele Diversity';
    plt.ylabel = '# Alleles / # Cells';
    dat.yvals = {apc.mu};
    plt.ylim = [0 0.6];

    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/AlleleDiversityOverTime', outdir));
    
    plt.title = 'Fraction of Cells Edited';
    plt.ylabel = 'Cells Edited (%)';    
    dat.yvals = {ee};
    plt.ylim = [0 100];
    
    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/FractionEditedOverTime', outdir));
    
    plt.title = 'CARLIN Potential';
    plt.ylabel = 'Unmodified Target Sites';
    dat.yvals = {CP.mu};
    plt.ylim = [0 10];
    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/CARLINPotentialOverTime', outdir));
    
    plt.title = 'Allele Diversity';    
    plt.ylabel = 'Diversity Index';    
    dat.yvals = {di.mu};
    plt.ylim = [0 0.5];
    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/CARLINDiversityIndexOverTime', outdir));
    
    close all;


end

