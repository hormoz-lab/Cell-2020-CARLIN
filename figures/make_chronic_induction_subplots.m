function make_chronic_induction_subplots(analysis_dir, results_dir)
   
    outdir = sprintf('%s/Figures/ChronicInduction', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    % Initial dataset for which we show edits/diversity and what is
    % mentioned in text
    
    exemplar = load(sprintf('%s/ChronicInduction/PosDox/Low/96h/Summary.mat', analysis_dir), 'summary');
    
    plot_indel_freq_vs_length(exemplar.summary);
    paper_print(sprintf('%s/LowDox96hInDelLengthFrequency', outdir));
    
    plot_highlighted_alleles(exemplar.summary, 50);
    paper_print(sprintf('%s/LowDox96hTop50Sequences', outdir));
    
    plot_site_decomposition(exemplar.summary, true, 'ES Cells', 'Cells');
    paper_print(sprintf('%s/LowDox96hSiteDecomposition', outdir));
        
    neg_folders = {'0h'; '24h'; '48h'; '72h'; '96h'};
    pos_folders = {'0h'; '12h'; '24h'; '48h'; '72h'; '96h'};
            
    N_neg = length(neg_folders);
    N_pos = length(pos_folders);
    
    chord = cell(4,1);
    chord{1} = load(sprintf('%s/ChronicInduction/PosDox/0h/Summary.mat', analysis_dir), 'summary');
    chord{2} = load(sprintf('%s/ChronicInduction/PosDox/Low/12h/Summary.mat', analysis_dir), 'summary');
    chord{3} = load(sprintf('%s/ChronicInduction/PosDox/Low/24h/Summary.mat', analysis_dir), 'summary');
    chord{4} = load(sprintf('%s/ChronicInduction/PosDox/Low/48h/Summary.mat', analysis_dir), 'summary');
        
    chord = cellfun(@(x) x.summary, chord, 'un', false);
    plot_quad_stargate(chord, pos_folders(1:4));
    paper_print(sprintf('%s/Stargate', outdir));
       
    chord_radius = 3;
    
    conc = {'Medium'; 'High'};
    tps = {'12h'; '24h'; '48h'};
    
    for i = 1:length(conc)
        for j = 1:length(tps)
            load(sprintf('%s/ChronicInduction/PosDox/%s/%s/Summary.mat', analysis_dir, conc{i}, tps{j}), 'summary');
             figure('Units', 'centimeters', 'Position', [0, 0, chord_radius, chord_radius], ...
                    'PaperUnits', 'centimeters', 'PaperSize', [chord_radius, chord_radius]);
             sp = plot_stargate.create(summary, 1, 1, 1);
             set(sp, 'Units', 'centimeters', 'Position', [0 0 chord_radius chord_radius]);
             axis tight;
             box off;
             paper_print(sprintf('%s/%s%sStargate', outdir, conc{i}, tps{j}));
        end
    end
    
    conc = {'Trial1'; 'Trial2'};
    tps = {'24h'; '48h'; '72h'};
    
    for j = 1:length(tps)
        for i = 1:length(conc)
            samples{i} = load(sprintf('%s/ChronicInduction/NegDox/%s/%s/Summary.mat', analysis_dir, conc{i}, tps{j}), 'summary');
        end
        summary = cellfun(@(x) x.summary, samples, 'un', false)';
        summary = ExperimentSummary.FromMerge(vertcat(summary{:}));
         figure('Units', 'centimeters', 'Position', [0, 0, chord_radius, chord_radius], ...
                'PaperUnits', 'centimeters', 'PaperSize', [chord_radius, chord_radius]);
         sp = plot_stargate.create(summary, 1, 1, 1);
         set(sp, 'Units', 'centimeters', 'Position', [0 0 chord_radius chord_radius]);
         axis tight;
         box off;
         paper_print(sprintf('%s/NegDox%sStargate', outdir, tps{j}));
    end
    
    ee = struct('Neg',    zeros(N_neg,1), ...
                'Low',    zeros(N_pos,1), ...
                'Medium', zeros(N_pos,1), ...
                'High',   zeros(N_pos,1));
    
    apc = struct('Neg',    struct('mu', zeros(N_neg,1), 'sig', zeros(N_neg,1)), ...
                 'Low',    struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                 'Medium', struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                 'High',   struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)));
        
    di = struct('Neg',    struct('mu', zeros(N_neg,1), 'sig', zeros(N_neg,1)), ...
                'Low',    struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                'Medium', struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                'High',   struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)));
    
    CP = struct('Neg',    struct('mu', zeros(N_neg,1), 'sig', zeros(N_neg,1), 'dist', {cell(N_neg,1)}), ...
                'Low',    struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1), 'dist', {cell(N_pos,1)}), ...
                'Medium', struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1), 'dist', {cell(N_pos,1)}), ...
                'High',   struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1), 'dist', {cell(N_pos,1)}));
    
    num_bp_del = struct('Neg',    struct('mu', zeros(N_neg,1), 'sig', zeros(N_neg,1)), ...
                        'Low',    struct('mu', zeros(N_pos,1), 'se', zeros(N_pos,1)), ...
                        'Medium', struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                        'High',   struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)));
                    
	num_bp_ins = struct('Neg',    struct('mu', zeros(N_neg,1), 'sig', zeros(N_neg,1)), ...
                        'Low',    struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                        'Medium', struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)), ...
                        'High',   struct('mu', zeros(N_pos,1), 'sig', zeros(N_pos,1)));

    edited_only = true;
  
    N_bootstrap_trials = 1000;
    rng(23849398);

    for i = 1:N_neg
        fprintf('Loading -Dox %d of %d\n', i, N_neg);
        con1 = load(sprintf('%s/ChronicInduction/NegDox/Trial1/%s/Summary.mat', analysis_dir, neg_folders{i}), 'summary');
        con2 = load(sprintf('%s/ChronicInduction/NegDox/Trial2/%s/Summary.mat', analysis_dir, neg_folders{i}), 'summary');        
        summary = ExperimentSummary.FromMerge([con1.summary; con2.summary]);
        ee.Neg(i) = sum(summary.allele_freqs(2:end))/sum(summary.allele_freqs)*100;
        [apc.Neg.mu(i), apc.Neg.sig(i)]               = bootstrap_alleles_per_cell(summary, N_bootstrap_trials);     
        [di.Neg.mu(i), di.Neg.sig(i)]                 = bootstrap_diversity_index(summary, N_bootstrap_trials);
        [CP.Neg.mu(i), CP.Neg.sig(i), CP.Neg.dist{i}] = CARLIN_potential(summary);
        [num_bp_del.Neg.mu(i), num_bp_ins.Neg.mu(i), num_bp_del.Neg.sig(i), num_bp_ins.Neg.sig(i)] = ...
            Mutation.num_bps_indel_stats(summary, edited_only);
    end
    
    conc = {'Low'; 'Medium'; 'High'};

    for i = 1:N_pos
        fprintf('Loading +Dox %d of %d\n', i, N_pos);
        if (i == 1)
            load(sprintf('%s/ChronicInduction/PosDox/%s/Summary.mat', analysis_dir, pos_folders{i}), 'summary');    
            
            ee.Low(i) = summary.N.eventful_tags/summary.N.called_tags*100;
            ee.Medium = ee.Low;
            ee.High = ee.Low;
            
            [apc.Low.mu(i), apc.Low.sig(i)] = bootstrap_alleles_per_cell(summary, N_bootstrap_trials);
            apc.Medium = apc.Low; 
            apc.High = apc.Low;
            
            [di.Low.mu(i), di.Low.sig(i)]  = bootstrap_diversity_index(summary, N_bootstrap_trials);
            di.Medium = di.Low;
            di.High = di.Low;
            
            [CP.Low.mu(i), CP.Low.sig(i), CP.Low.dist{i}] = CARLIN_potential(summary);
            CP.Medium = CP.Low;
            CP.High = CP.Low;
            
            [num_bp_del.Low.mu(i), num_bp_ins.Low.mu(i), num_bp_del.Low.sig(i), num_bp_ins.Low.sig(i)] = ...
                Mutation.num_bps_indel_stats(summary, edited_only);
            num_bp_del.Medium = num_bp_del.Low;
            num_bp_del.High   = num_bp_del.Low;
            num_bp_ins.Medium = num_bp_ins.Low;
            num_bp_ins.High   = num_bp_ins.Low;
            
        else
            for j = 1:length(conc)                
                load(sprintf('%s/ChronicInduction/PosDox/%s/%s/Summary.mat', analysis_dir, conc{j}, pos_folders{i}), 'summary');
                ee.(conc{j})(i) = summary.N.eventful_tags/summary.N.called_tags*100;
                [apc.(conc{j}).mu(i), apc.(conc{j}).sig(i)] = bootstrap_alleles_per_cell(summary, N_bootstrap_trials);
                [di.(conc{j}).mu(i) , di.(conc{j}).sig(i) ] = bootstrap_diversity_index(summary, N_bootstrap_trials);
                [CP.(conc{j}).mu(i), CP.(conc{j}).sig(i), CP.(conc{j}).dist{i}] = CARLIN_potential(summary);
                [num_bp_del.(conc{j}).mu(i), num_bp_ins.(conc{j}).mu(i), num_bp_del.(conc{j}).sig(i), num_bp_ins.(conc{j}).sig(i)] = ...
                    Mutation.num_bps_indel_stats(summary, edited_only);
            end
        end
    end
        
    dat.series_names = {'No Dox'; 'Low Dox'; 'Medium Dox'; 'High Dox'};
    dat.xvals = [{[1; 3; 4; 5; 6]}; repmat({[1:6]'}, 3, 1)];    
    dat.markers = {'^'; 's'; 'o'; 'v'};
    
    plt.xlim = [0.5 6.5];    
    plt.xticklabels = pos_folders;    
    plt.xlabel = 'Time (h)';
    
    plt.title = 'BPs Deleted In Edited Cells';
    plt.ylabel = 'Average BPs Deleted';    
    dat.yvals = struct2cell(structfun(@(x) x.mu, num_bp_del, 'un', false));
    plt.ylim = [0 250];

    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/BPsDeletedOverTime', outdir));
    
    plt.title = 'Allele Diversity';
    plt.ylabel = '# Alleles / # Cells';
    dat.yvals = struct2cell(structfun(@(x) x.mu, apc, 'un', false));
    plt.ylim = [0 0.6];

    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/AlleleDiversityOverTime', outdir));
    
    plt.title = 'Fraction of Cells Edited';    
    plt.ylabel = 'Cells Edited (%)';    
    dat.yvals = struct2cell(ee);
    plt.ylim = [0 100];
    
    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/FractionEditedOverTime', outdir));
    
    plt.title = 'CARLIN Potential';
    plt.ylabel = 'Unmodified Target Sites';
    dat.yvals = struct2cell(structfun(@(x) x.mu, CP, 'un', false));
    plt.ylim = [0 10];
    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/CARLINPotentialOverTime', outdir));
    
    plt.title = 'Allele Diversity';    
    plt.ylabel = 'Diversity Index';    
    dat.yvals = struct2cell(structfun(@(x) x.mu, di, 'un', false));
    plt.ylim = [0 0.5];
    plot_errorbar_multiseries(dat, plt);
    paper_print(sprintf('%s/CARLINDiversityIndexOverTime', outdir));
    
    plot_errorbar_multiseries(dat, plt, true);
    paper_print(sprintf('%s/Legend', outdir));
    
    close all;

end
