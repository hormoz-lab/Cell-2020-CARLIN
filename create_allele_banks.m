function create_allele_banks(analysis_path, results_path)
    
    fprintf('Creating Protocol 1 +Dox Bank...\n');
    
    FO215 = combine_samples(analysis_path, 'FO215', {'1'; '2'; '3'});
    FO217 = combine_samples(analysis_path, 'FO217', {'1'; '2'; '3'});
    FO395 = combine_samples(analysis_path, 'FO395', {'125K'; '25K'; '5K'; '1K'});
    FO397 = combine_samples(analysis_path, 'FO397', {'125K'; '25K'; '5K'; '1K'});
    
    Bank.Create([FO215; FO217; FO395; FO397], {'FO215'; 'FO217'; 'FO395'; 'FO397'}, sprintf('%s/Banks/Protocol1/RNA/PosDox', results_path));
    
    fprintf('Creating Protocol 1 -Dox Bank...\n');
    
    FOND1 = combine_samples(analysis_path, 'FOND1', {'125K'; '25K'; '5K'; '1K'});
    FOND2 = combine_samples(analysis_path, 'FOND2', {'125K'; '25K'; '5K'; '1K'});
    
    Bank.Create([FOND1; FOND2], {'FOND1'; 'FOND2'}, sprintf('%s/Banks/Protocol1/RNA/NegDox', results_path));
        
    fprintf('Creating Protocol 1 -Cas9 Bank...\n');    
  
    fprintf('...FO803\n');
    dat = load(sprintf('%s/Diversity/FO803/Summary.mat', analysis_path), 'summary');
    Bank.Create(dat.summary, {'FO803'}, sprintf('%s/Banks/Protocol1/RNA/NegCas9', results_path));
    clear dat;
    
    fprintf('Creating Protocol 2 +Dox Bank...\n');
    
    FO912 = combine_samples(analysis_path, 'FO912', {'L/125K'; 'L/25K'; 'L/5K'; 'R/125K'; 'R/25K'; 'R/5K'});
    FO913 = combine_samples(analysis_path, 'FO913', {'L/125K'; 'L/25K'; 'L/5K'; 'R/125K'; 'R/25K'; 'R/5K'});
    FO914 = combine_samples(analysis_path, 'FO914', {'L/125K'; 'L/25K'; 'L/5K'; 'R/125K'; 'R/25K'; 'R/5K'});
       
    Bank.Create([FO912; FO913; FO914], {'FO912'; 'FO913'; 'FO914'}, sprintf('%s/Banks/Protocol2/RNA/PosDox', results_path));
    
    fprintf('Creating Protocol 2 -Dox Bank...\n');
    
    FO915 = combine_samples(analysis_path, 'FO915', {'L/125K'; 'L/25K'; 'L/5K'; 'R/125K'; 'R/25K'; 'R/5K'});
    FO916 = combine_samples(analysis_path, 'FO916', {'L/125K'; 'L/25K'; 'L/5K'; 'R/125K'; 'R/25K'; 'R/5K'});    
    
    Bank.Create([FO915; FO916], {'FO915'; 'FO916'}, sprintf('%s/Banks/Protocol2/RNA/NegDox', results_path));
    
    fprintf('Pooling bulk samples for test replicates...\n');
    
    SB213 = combine_samples(analysis_path, 'SB213', {'1'; '2'; '3'; '4'});
    SB214 = combine_samples(analysis_path, 'SB214', {'1'; '2'; '3'; '4'});
   
end 

function summary = combine_samples(analysis_path, sample_path, split_path)
    
    fprintf('...%s\n', sample_path);

    sample_names = cellfun(@(x) sprintf('%s/%s', sample_path, x), split_path, 'un', false);
    dat = cell(length(sample_names),1);
    for i = 1:length(sample_names)
        dat{i} = load(sprintf('%s/Diversity/%s/Summary.mat', analysis_path, sample_names{i}), 'summary');
    end
    dat = vertcat(dat{:});
    [summary, sample_map, allele_breakdown_by_sample] = ExperimentSummary.FromMerge(vertcat(dat(:).summary));
    outpath = sprintf('%s/Diversity/%s/Combined', analysis_path, sample_path);
    if (~exist(outpath, 'dir'))
        mkdir(outpath);
    end
    save(sprintf('%s/Summary.mat', outpath), 'summary', 'sample_map', 'allele_breakdown_by_sample', 'sample_names');
    plot_summary(summary, outpath);
end
