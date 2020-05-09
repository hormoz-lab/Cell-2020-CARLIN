function run_paper_experiments(experiment_series, processed_path, analysis_path, cmd_path)

    if (nargin < 4)
        evaluate = true;
    else
        evaluate = false;        
        if (~exist(cmd_path, 'dir'))
            mkdir(cmd_path);
        end
        fid = fopen(sprintf('%s/%s.txt', cmd_path, experiment_series), 'wt');
    end
    
    if (strcmp(experiment_series, 'Sanger'))
            
        cfg = parse_config_file('Sanger');
        [folder, ~, ~] = fileparts(mfilename('fullpath'));
        orig_clones = upper(splitlines(fileread(sprintf('%s/data/Sanger/Sanger8.txt', folder))));
        [trimmed_clones, reorder] = FastQData.extract_CARLIN_from_sequences(orig_clones, [1:size(orig_clones,1)]', cfg);
        [~, trimmed_clones] = cellfun(@(x) CARLIN_def.cas9_align(x), trimmed_clones(reorder), 'un', false);
        trimmed_clones = AlignedSEQ.sanitize_prefix_postfix(trimmed_clones);
        trimmed_clones = AlignedSEQ.sanitize_conserved_regions(trimmed_clones);
        for i = 1:size(trimmed_clones,1)
            summary = ExperimentSummary(trimmed_clones(i));
            savedir = sprintf('%s/InVitroPhylogeny/Gen1/%d', analysis_path, i);
            if (~exist(savedir, 'dir'))
                mkdir(savedir);        
            end
            save(sprintf('%s/Summary.mat', savedir), 'summary');
            plot_summary(summary, savedir);
            Mutation.ToFile(summary.alleles, savedir, 'AlleleAnnotations.txt');
        end  
        
    else
        
        [sample_sheet, samples_to_run] = get_experiment_metadata(experiment_series);
    
        if (strcmp(experiment_series, 'Transcriptome'))

            fprintf('Preparing transcriptome for %d single cell experiments\n', length(samples_to_run));
            for s = samples_to_run
                [input_file, output_path] = prep_transcriptome_inputs(sample_sheet, s, processed_path, analysis_path);
                cmd = sprintf('process_transcriptome(''%s'', ''%s'');', input_file, output_path);
                if (evaluate)
                    eval(cmd);
                else
                    fprintf(fid, '%s\n', strrep(cmd, '''', '\'''));
                end
            end

        elseif (ismember(experiment_series, {'5FU'; 'EB'; 'SCReplicate'}))

            fprintf('Running CARLIN amplicon pipeline on %d samples for %s experiment with corresponding transcriptome data\n', ...
                length(samples_to_run), experiment_series);
            for s = samples_to_run
                [input_file, cfg_type, output_path, CB_file] = prep_SC_inputs(sample_sheet, s, processed_path, analysis_path);                
                cmd = sprintf('analyze_CARLIN(%s, ''%s'', ''%s'', ''ref_CB_file'', ''%s'', ''read_override_UMI_denoised'', 10, ''read_override_CB_denoised'', 10);', ...
                              stringify(input_file), cfg_type, output_path, CB_file);
                if (evaluate)
                    assert(exist(CB_file, 'file')==2, 'Process transcriptome before running %s', experiment_series);
                    eval(cmd);
                else
                    fprintf(fid, '%s\n', strrep(cmd, '''', '\'''));
                end
            end

        else

            fprintf('Running CARLIN amplicon pipeline on %d samples for %s experiment\n', length(samples_to_run), experiment_series);    
            for s = samples_to_run
                [input_file, cfg_type, output_path, num_cells] = prep_bulk_inputs(sample_sheet, s, processed_path, analysis_path);
                if (num_cells == '?')
                    cmd = sprintf('analyze_CARLIN(''%s'', ''%s'', ''%s'');', input_file, cfg_type, output_path);                    
                else
                    cmd = sprintf('analyze_CARLIN(''%s'', ''%s'', ''%s'', ''max_molecules'', %d);', input_file, cfg_type, output_path, num_cells);                    
                end
                if (evaluate)
                    eval(cmd);
                else
                    fprintf(fid, '%s\n', strrep(cmd, '''', '\'''));
                end
            end 
        end
    end
    
    if (~evaluate)
        fclose(fid);
    end
end

function s = stringify(files)
    s = '{';
    for i = 1:size(files,1)
        for j = 1:size(files,2)
            if (j > 1)
                s = [s ', '];
            end
            s = [s '''' files{i,j} ''''];
        end
        if (i < size(files,1))
            s = [s '; '];
        end
    end
    s = [s '}'];
end