function make_5FU_subplots(results_dir)

    datfile = sprintf('%s/5FU/Analysis.mat', results_dir);
    assert(exist(datfile, 'file')==2);
    load(datfile);
    
    outdir = sprintf('%s/Figures/5FU', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    title_str = {'5-FU Mouse 2'; '5-FU Mouse 1'; 'Control 2'; 'Control 1'; 'Control 3'};
    
    ceres_names = {'Parent HSC'; 'Childless HSC'; 'HSC-Derived Progeny'; 'Progeny-Only Clone'};    
    
    color = [0 0.7 0;
             0 0 1;
             0 0.7 0;
             0.7 0.7 0];
          
    marker = ['oodo']';
    
    rng(304394839);
    for i = [1 3 4 5]
        replicate_equivalence(data_path{i}, title_str{i}, outdir, i==1);
    end
    
    plot_dotplot(dp_data, marker_genes, phenotypes, 'Marker Gene Expression in Clusters');
    paper_print(sprintf('%s/DotPlot', outdir));

    highlight_patch = get_highlight_patch(pooled.x, pooled.y, find(pooled.cg_pheno==hsc_pheno));    
    highlight_patch = [pooled.x(highlight_patch) pooled.y(highlight_patch)];
    
    for i = [pos_samples(1:2) neg_samples(1:2)]
        
        ceres = NaN(size(samples{i}.x));
        ceres(samples{i}.cg_pheno == hsc_pheno & ismember(samples{i}.allele, samples{i}.sig_alleles.hsc_derived_allele )) = 1;
        ceres(samples{i}.cg_pheno == hsc_pheno & ismember(samples{i}.allele, samples{i}.sig_alleles.hsc_only_allele    )) = 2;
        ceres(samples{i}.cg_pheno ~= hsc_pheno & ismember(samples{i}.allele, samples{i}.sig_alleles.hsc_derived_allele )) = 3;
        ceres(samples{i}.cg_pheno ~= hsc_pheno & ismember(samples{i}.allele, samples{i}.sig_alleles.derived_only_allele)) = 4;
        
        temp_title = sprintf('%s - %d sig. clones with %d cells', title_str{i}, ...
            length(samples{i}.sig_alleles.all), sum(samples{i}.summary.allele_freqs(samples{i}.sig_alleles.all)));
        
        plot_allele_on_transcriptome('5FU', samples{i}.x, samples{i}.y, ceres, ceres_names, [12; 12; 12; 3], ...
                                     color, marker, [4; 3; 2; 1], true, [1; 3; 2; 4], temp_title, false, highlight_patch);
        paper_print(sprintf('%s/ClonalExpansion%sWithLegend', outdir, mouse{i}));
        
        plot_allele_on_transcriptome('5FU', samples{i}.x, samples{i}.y, ceres, ceres_names, [12; 12; 12; 3], ...
                                     color, marker, [4; 3; 2; 1], false, [1; 3; 2; 4], temp_title, false, highlight_patch);
        paper_print(sprintf('%s/ClonalExpansion%sWithoutLegend', outdir, mouse{i}));
        
    end    
  
    plot_phenotype_on_transcriptome(pooled.x, pooled.y, pooled.cg_pheno, coarse_grain_pheno(:,1), ...
                                    [3; 2; 5; 4; 1; 6; 7], 'Cluster-Based Cell Annotations', true, false);
    paper_print(sprintf('%s/Phenotypes', outdir));
    
    plot_order = cellfun(@fliplr, phenotypes(:,2), 'un', false);
    plot_order = fliplr(horzcat(plot_order{:})+1);
    
    plot_phenotype_on_transcriptome(pooled.x, pooled.y, pooled.louvain+1, cellstr(num2str([0:max(pooled.louvain)]')), ...
                                    plot_order', 'Louvain Clusters', false, true);
    paper_print(sprintf('%s/LouvainClusters', outdir));
    
    plot_phenotype_on_transcriptome(pooled.x, pooled.y, pooled.sample, title_str, ...
                                    [4; 1; 5; 2; 3], 'Aligned Samples', false, false, [2 1 4 3 5]);
    paper_print(sprintf('%s/SampleSplit', outdir));
    
    load(sprintf('%s/SC/Stats.mat', results_dir), 'FU5');
    
    n_cells = cellfun(@(x) sum(ismember(x.allele, union(x.sig_alleles.hsc_only_allele, x.sig_alleles.hsc_derived_allele))), samples);
    
    for i = 1:2
    
        pstring = get_pstring(FU5.acs{i,i}.pval);    

        xlabels = {sprintf('%s?  (n=%d)  ',      title_str{neg_samples(i)}, n_cells(neg_samples(i))); 
                   sprintf(' %s ?       (n=%d)', title_str{pos_samples(i)}, n_cells(pos_samples(i)))};
        xlabels = cellfun(@(x) strrep(x,'?','\newline'), xlabels,'UniformOutput',false);

        plot_twin_violin(nonzeros(FU5.acs{i,i}.neg.clone_size_freqs), find(FU5.acs{i,i}.neg.clone_size_freqs), ...
                         nonzeros(FU5.acs{i,i}.pos.clone_size_freqs), find(FU5.acs{i,i}.pos.clone_size_freqs), xlabels, pstring);
        paper_print(sprintf('%s/SCViolin_%s_%s', outdir, mouse{pos_samples(i)}, mouse{neg_samples(i)}));
        
        plot_number_clones_downsampled(samples([neg_samples(i), pos_samples(i)]), color([1 2 4],:), ...
                                       title_str([neg_samples(i), pos_samples(i)]), '# Clones');
        paper_print(sprintf('%s/NumberClones_%s_%s', outdir, mouse{pos_samples(i)}, mouse{neg_samples(i)}));
        
        if (i == 2)
            plot_5FU_clone_pheno_heatmap(samples([neg_samples(i) pos_samples(i)]), coarse_grain_pheno(:,1), ...
                                         title_str([neg_samples(i) pos_samples(i)]), color([1,2,4],:));
            paper_print(sprintf('%s/Heatmap_%s_%s', outdir, mouse{pos_samples(i)}, mouse{neg_samples(i)}));
        else
            plot_5FU_clone_pheno_heatmap(samples([neg_samples(i) pos_samples(i)]), coarse_grain_pheno(:,1), ...
                                         title_str([neg_samples(i) pos_samples(i)]), color([1,2,4],:), length(samples{pos_samples(2)}.sig_alleles.all));
            paper_print(sprintf('%s/Heatmap_%s_%s', outdir, mouse{pos_samples(i)}, mouse{neg_samples(i)}));
        end
    end
    
    close all;

end

function replicate_equivalence(data_path, title_str, outdir, legend)

    N_trials = 10000;
    CI = 0.95;

    load(sprintf('%s/Amplicon/Replicate/Combined/Summary.mat', data_path), 'combined');    
    is_template = strcmp(cellfun(@(x) degap(x.get_seq), combined.summary.alleles, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
    [~, KS_dist] = bootstrap_KS_statistic(combined.allele_breakdown_by_sample, N_trials);
    KS_report = percentile(KS_dist, CI);
    
    plot_replicate_proportions(combined.allele_breakdown_by_sample, is_template, sprintf('%s: KS_{%.2f}=%.1e', title_str, CI, KS_report), legend);
    paper_print(sprintf('%s/Replicate%s', outdir, deblank(title_str)));
    
end

function pstring = get_pstring(pval)
    if (pval < 1e-6)
        pstring = 'p < 1e-6';
    elseif (pval >= 0.1)
        pstring = sprintf('p = %.1g', pval);
    else
        pstring = sprintf('p = %2.1e', pval);
    end
end
