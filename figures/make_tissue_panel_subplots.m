function make_tissue_panel_subplots(analysis_dir, results_dir)
    
    outdir = sprintf('%s/Figures/TissuePanel', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    samples = {'FO912'; 'FO916'};
    tissue_labels = {'Skin'; 'Brain'; 'Granulocytes'; 'Muscle'; 'Heart'; 'Lung'; 'Liver'; 'Intestine'; 'Spleen'; 'Kidney'; 'Gonad'};
    N_samples = length(samples);
    N_tissues = length(tissue_labels);
    ee = zeros(N_tissues, N_samples);
    
    for i = 1:N_samples
        for j = 1:N_tissues            
            if (~strcmp(tissue_labels{j}, 'Granulocytes'))
                load(sprintf('%s/TissuePanel/%s/%s/Summary.mat', analysis_dir, samples{i}, tissue_labels{j}), 'summary');
                ee(j,i) = summary.N.eventful_tags / summary.N.called_tags * 100;
                if (i == 1)
                    plot_site_decomposition(summary, false, tissue_labels{j}, 'Transcripts');
                end
            else
                if (i == 1)
                    load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir));
                    summary = bank.summary;
                else
                    load(sprintf('%s/Banks/Protocol2/RNA/NegDox/Bank.mat', results_dir));
                    summary = bank.summary;
                end                
                is_template = strcmp(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
                ee(j,i) = sum(summary.allele_freqs(~is_template)) / sum(summary.allele_freqs) * 100;
                if (i == 1)
                    plot_site_decomposition(summary, true, tissue_labels{j}, 'Transcripts');
                end
            end
            if (i == 1)
                paper_print(sprintf('%s/%s', outdir, tissue_labels{j}));            
            end
        end
    end
    
    plot_tissue_panel_summary(tissue_labels, ee, 'Fraction of Transcripts Edited', 'Transcripts edited (%)', {'+Dox'; '-Dox'}, true);
    paper_print(sprintf('%s/TissuePanelFractionEdited', outdir));
    

    tissue_labels = {'Lung'; 'BM'; 'Skin'; 'Liver'; 'Intestine'};    
    cell_labels = {{'Monocytes'; 'Epithelial'; 'Mesenchymal'}; 
                   {'Monocytes'; 'Mesenchymal'}; 
                   {'Epithelial'};
                   {'Mesenchymal'};
                   {'Epithelial'}};
    
    N_tissues = length(tissue_labels);
    ee = cell(N_tissues, 1);
    
    mouse_set1 = {'SB213'; 'SB214'};
    mouse_set2 = {'SB225'; 'SB226'};
    
    for i = 1:N_tissues
        if (strcmp(tissue_labels{i}, 'Liver'))
            mouse_set = mouse_set2;
        else
            mouse_set = mouse_set1;
        end
        N_mice = length(mouse_set);
        ee{i} = zeros(length(cell_labels{i}), N_mice);
        for j = 1:length(cell_labels{i})
            for k = 1:N_mice
                load(sprintf('%s/TissuePanel/%s/%s/%s/Summary.mat', analysis_dir, mouse_set{k}, tissue_labels{i}, cell_labels{i}{j}), 'summary');
                ee{i}(j,k) = summary.N.eventful_tags / summary.N.called_tags * 100;
            end
        end
    end
    
    plot_tissue_panel_summary(cell_labels, ee, '               Lung                           BM            Skin    Liver  Intestine', ...
                              'Transcripts edited (%)', {'Replicate 1'; 'Replicate 2'}, false);    
    paper_print(sprintf('%s/CellPanelFractionEdited', outdir));

    close all;
    
end