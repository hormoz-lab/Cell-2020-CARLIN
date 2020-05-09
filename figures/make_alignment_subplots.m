function make_alignment_subplots(analysis_dir, results_dir)

    outdir = sprintf('%s/Figures/Alignment', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    ref = CARLIN_def.getInstance;
    L = ref.width.CARLIN;
    trim_amount = ref.width.Primer5+ref.width.Primer3+ref.width.SecondarySequence;
    min_length = ref.width.min_length-0.5;
    max_length = 300;
    
    % Count anything with 1.5 bps of template as uncut. Fragment analysis
    % can't seem to resolve finer than this based on loadings even with
    % -Cas9 samples
    envelope = [1.5 1.5];
    
    [metadata, samples_to_run] = get_experiment_metadata('ChronicInduction');
    metadata = metadata(samples_to_run,:);
   
    % Low Dox @72h fragment analysis is corrupted, and looks bad aesthetically
    % to include the time point for only 2/3 concentrations so just ignore
    % it altogether. 0h data is common for all three concentrations so
    % makes it annoying to plot on grid. Just ignore 
    metadata = metadata(contains(metadata.ProcessedPath, 'PosDox') & ~contains(metadata.ProcessedPath, {'72h'; '0h'}),:);
    
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    fragment_file = sprintf('%s/../data/FragmentAnalysis/ChronicInduction.xlsx', folder);
    
    for i = 1:size(metadata,1)
        T = parse_fragment_analysis_spreadsheet(fragment_file, sprintf('O%s', metadata.SampleIndex{i}));
        load(sprintf('%s/%s/Summary.mat', analysis_dir, metadata.ProcessedPath{i}), 'summary');    
        [frag.xvals, frag.CDF, frag.fraction_edited] = process_fragment_analysis(T, trim_amount, min_length, max_length, L, envelope);
        [algo.xvals, algo.CDF, algo.fraction_edited] = length_distribution_from_summary(summary, min_length, max_length);
        
        if (contains(metadata.ProcessedPath{i}, 'Low') && contains(metadata.ProcessedPath{i}, '12h'))
            show_legend = true;
        else
            show_legend = false;
        end
        
        title_str = strrep(erase(metadata.ProcessedPath{i}, 'ChronicInduction/PosDox/'), '/', ' Dox @ ');

        plot_fragment_analysis_CDF(frag, algo, L, envelope, title_str, show_legend);
        paper_print(sprintf('%s/%s', outdir, deblank(title_str)));
    end
    
    close all;
    
end

function [xvals, CDF, fraction_edited] = length_distribution_from_summary(summary, min_length, max_length)

    [mid_bp, ~, freqs] = find(accumarray(cellfun(@(x) length(degap(x.get_seq)), summary.alleles), summary.allele_freqs));
    freqs = cumsum(freqs);
    
    xvals = reshape([mid_bp-0.5 mid_bp+0.5]', [2*numel(mid_bp), 1]);
    CDF = zeros(size(xvals));    
    CDF(2:2:end) = freqs;
    xvals = [min_length; xvals];
    CDF = [0; CDF];
    CDF(2:2:end) = CDF(1:2:end-1);
    last_ind = find(xvals <= max_length, 1, 'last');
    xvals = [xvals(1:last_ind,:); max_length];
    CDF   = [CDF(1:last_ind); CDF(end)];
    CDF   = CDF/CDF(end);
    assert(length(xvals) == length(CDF));
    assert(issorted(xvals) && issorted(CDF));
    
    is_template = strcmp(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
    fraction_edited = sum(summary.allele_freqs(~is_template))/sum(summary.allele_freqs);
    
end
