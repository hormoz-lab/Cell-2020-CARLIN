function pool_samples(analysis_dir)

    fprintf('Pooling data for FO864\n');    
    samples = cell(4,1);
    samples{1} = load(sprintf('%s/EB/FO864/SC/LF/Amplicon/Summary.mat', analysis_dir));
    samples{2} = load(sprintf('%s/EB/FO864/SC/RF/Amplicon/Summary.mat', analysis_dir));
    samples{3} = load(sprintf('%s/EB/FO864/SC/LH/Amplicon/Summary.mat', analysis_dir));
    samples{4} = load(sprintf('%s/EB/FO864/SC/RH/Amplicon/Summary.mat', analysis_dir));   
    sample_names = {'FO864/LF'; 'FO864/RF'; 'FO864/LH'; 'FO864/RH'};
    summary = cellfun(@(x) ExperimentSummary(x.summary), samples, 'un', false);
    [combined.summary, combined.sample_map, combined.allele_breakdown_by_sample] = ExperimentSummary.FromMerge(vertcat(summary{:}));
    outdir = sprintf('%s/EB/FO864/SC/Combined/Amplicon', analysis_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    save(sprintf('%s/Summary.mat', outdir), 'combined', 'sample_names', 'samples');
    plot_summary(combined.summary, outdir);
    
    fprintf('Pooling data for 5FU Replicates\n');
    sample_names = {'Rep1'; 'Rep2'};
    data_paths = {'FO898/SC'; 'FO892/SC/1'; 'FO892/SC/2'; 'FO897/SC'};
    for i = 1:length(data_paths)
        samples = cell(2,1);     
        samples{1} = load(sprintf('%s/5FU/%s/Amplicon/Replicate/1/Summary.mat', analysis_dir, data_paths{i}));
        samples{2} = load(sprintf('%s/5FU/%s/Amplicon/Replicate/2/Summary.mat', analysis_dir, data_paths{i}));        
        summary = cellfun(@(x) ExperimentSummary(x.summary), samples, 'un', false);
        [combined.summary, combined.sample_map, combined.allele_breakdown_by_sample] = ExperimentSummary.FromMerge(vertcat(summary{:}));
        outdir = sprintf('%s/5FU/%s/Amplicon/Replicate/Combined', analysis_dir, data_paths{i});
        if (~exist(outdir, 'dir'))
            mkdir(outdir);
        end
        save(sprintf('%s/Summary.mat', outdir), 'combined', 'sample_names', 'samples');
        plot_summary(combined.summary, outdir);
    end

end