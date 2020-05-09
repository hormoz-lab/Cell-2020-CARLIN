function make_inducible_subplots(analysis_dir, results_dir)
    
    outdir = sprintf('%s/Figures/Inducible', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    [sample_sheet, samples_to_run] = get_experiment_metadata('Inducible');
    data_paths = sample_sheet.ProcessedPath(samples_to_run);
    N = length(data_paths);
    summary = cell(N, 1);
    
    for i = 1:N
        summary{i} = load(sprintf('%s/%s/Summary.mat', analysis_dir, data_paths{i}), 'summary');
    end
    summary = cellfun(@(x) x.summary, summary, 'un', false);
    
    FO409 = ExperimentSummary.FromMerge(vertcat(summary{1:4}));
    FO406 = ExperimentSummary.FromMerge(vertcat(summary{5:8}));    
    [summary, ~, allele_breakdown_by_sample] = ExperimentSummary.FromMerge([FO409; FO406]);
    
    plot_highlighted_alleles(summary, 50);
    paper_print(sprintf('%s/Top50Sequences', outdir));
    
    plot_allele_frequency_CDF(summary, 'iCARLIN');
    paper_print(sprintf('%s/CDF', outdir));
    
    plot_site_decomposition(summary, false, '', 'Transcripts');
    paper_print(sprintf('%s/SiteDecomposition', outdir));
    
    plot_indel_freq_vs_length(summary, allele_breakdown_by_sample);
    paper_print(sprintf('%s/InDelFreqVsLength', outdir));
    
    close all;

end