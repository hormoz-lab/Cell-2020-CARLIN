function generate_supplemental_tables(processed_path, analysis_path, results_path)

    outdir = sprintf('%s/SupplementalTables', results_path);
    
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
%     Supplementary Table 3
%     A reviewer asked for the number of reads corrected, which is
%     unfortunately not recorded in the Summary struct, so I have to unzip
%     and load the Analysis.mat files which is slow. 
    
    [sample_sheet, samples_to_run] = get_experiment_metadata('Diversity');
    samples = {'FO912'; 'FO913'; 'FO914'};
    read_summary = cell(length(samples),1);
    for s = 1:length(samples)
        samples_to_run_pos = samples_to_run(contains(sample_sheet.ProcessedPath(samples_to_run), samples{s}));
        summary = cell(length(samples_to_run_pos),1);
        corrected = zeros(length(samples_to_run_pos),1);
        for i = 1:length(samples_to_run_pos)
            [~, ~, output_path, ~] = prep_bulk_inputs(sample_sheet, samples_to_run_pos(i), processed_path, analysis_path);
            summary{i} = load(sprintf('%s/Summary.mat', output_path), 'summary');
            unzip_path = tempname;
            fprintf('Unzipping Analysis file to: %s/%s\n', unzip_path, 'Analysis.mat');
            gunzip(sprintf('%s/Analysis.mat.gz', output_path), unzip_path);
            assert(exist(sprintf('%s/Analysis.mat', unzip_path), 'file') == 2, sprintf('Unzip unsuccessful'));
            load(sprintf('%s/Analysis.mat', unzip_path), 'aligned', 'FQ');        
            delete(sprintf('%s/Analysis.mat', unzip_path));
            corrected(i) = dot(cellfun(@(x,y) ~isequal(degap(x.get_seq), y), aligned.aligned_SEQ(aligned.alignment_map), aligned.unaligned_SEQ), ...
                               accumarray(FQ.read_SEQ_valid(FQ.masks.valid_lines), 1));

        end
        summary = cellfun(@(x) x.summary.reads, summary, 'un', false);
        summary = vertcat(summary{:});        
        read_summary{s} = [cellfun(@(x) sum(vertcat(summary.(x))), fieldnames(summary)); sum(corrected)];        
    end
    
    T = table(read_summary{1}, read_summary{2}, read_summary{3}, 'VariableNames', samples, 'RowNames', [fieldnames(summary); 'corrected']);
    writetable(T, sprintf('%s/BankSummary.txt', outdir), 'Delimiter', '\t', 'WriteRowNames', true);
    
    % Supplementary Table 4
    
    paths = {'5FU/FO892/SC/2';
             '5FU/FO892/SC/1';
             '5FU/FO897/SC';
             '5FU/SB133/SC';
             '5FU/FO898/SC';             
             '5FU/FO817/SC';
             '5FU/FO837/SC';
             'EB/FO864/SC/LF';
             'EB/FO864/SC/RF';
             'EB/FO864/SC/LH';
             'EB/FO864/SC/RH'};
         
    quants = {'Sample';
              'MinimumUMICount'; 
              'FilteredTranscriptomeCBs';
              'MedianGenes';
              'FilteredAmpliconCBs';
              'MeanReads';
              'Alleles';
              'PctEdited'};
          
    T = table('Size', [length(paths), length(quants)], 'VariableTypes', [{'string'}; repmat("double", length(quants)-1, 1)], ...
              'VariableNames', quants, 'RowNames', paths);
    for i = 1:length(paths)
        load(sprintf('%s/%s/Transcriptome/Transcriptome.mat', analysis_path, paths{i}));
        T.Sample(i) = paths{i};
        T.MinimumUMICount(i) = ceil(transcriptome.cutoffs.umi);
        T.MedianGenes(i) = median(sum(transcriptome.counts(transcriptome.masks.cells.joint,:)>0,2));
        load(sprintf('%s/%s/Amplicon/Summary.mat', analysis_path, paths{i}));
        T.FilteredTranscriptomeCBs(i) = length(ref_CBs);        
        T.FilteredAmpliconCBs(i) = summary.N.called_tags;
        T.MeanReads(i) = round(summary.reads.called_tags_total/summary.N.called_tags);
        T.Alleles(i) = length(summary.alleles);
        T.PctEdited(i) = 100*summary.N.eventful_tags/summary.N.called_tags;
    end
        
    writetable(T, sprintf('%s/SCSummary.txt', outdir), 'Delimiter', '\t');
    
    % Supplementary Table 5

    load(sprintf('%s/SC/Stats.mat', results_path), 'FU5', 'EB');
    
    dlmwrite(sprintf('%s/NumberClones.txt', outdir), cellfun(@(x) x.pval, [EB.nc; FU5.nc]));
    dlmwrite(sprintf('%s/AverageHSCCloneSize.txt', outdir), cellfun(@(x) x.pval, [EB.acs; FU5.acs]));
    dlmwrite(sprintf('%s/NonUniformHSCCloneSize.txt', outdir), cellfun(@(x) x.pval, [{EB.csd}, FU5.csd]));    
    
end