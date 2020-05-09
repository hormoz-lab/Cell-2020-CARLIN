function make_EB_subplots(results_dir)

    datfile = sprintf('%s/EB/Analysis.mat', results_dir);
    assert(exist(datfile, 'file')==2);
    load(datfile);
    
    outdir = sprintf('%s/Figures/EB', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    plot_dotplot(dp_data, marker_genes, phenotypes, 'Marker Gene Expression in Clusters');
    paper_print(sprintf('%s/Figures/EB/DotPlot', results_dir));
    
    plot_phenotype_on_transcriptome(pooled.y, -pooled.x, pooled.cg_pheno, coarse_grain_pheno(:,1), ...
                                    [3; 2; 5; 4; 1; 6; 7], 'Cluster-Based Cell Annotations', true, false);
    paper_print(sprintf('%s/Figures/EB/Phenotypes', results_dir));
    
    plot_order = cellfun(@fliplr, phenotypes(:,2), 'un', false);
    plot_order = fliplr(horzcat(plot_order{:})+1);
    
    plot_phenotype_on_transcriptome(pooled.y, -pooled.x, pooled.louvain+1, cellstr(num2str([0:max(pooled.louvain)]')), ...
                                    plot_order', 'Louvain Clusters', false, true);
    paper_print(sprintf('%s/Figures/EB/LouvainClusters', results_dir));
    
    title_str = {'Left Femur'; 'Right Femur'; 'Left Humerus'; 'Right Humerus'};
    plot_phenotype_on_transcriptome(pooled.y, -pooled.x, pooled.sample, title_str, ...
                                    [4; 1; 2; 3], 'Aligned Bones', false, false);    
    paper_print(sprintf('%s/Figures/EB/SampleSplit', results_dir));
   
    pstring1 = get_pstring(allele_equal_rep, bb_pval(allele_equal_rep)*bonferroni_correction_factor);
    pstring2 = get_pstring(allele_biased_rep, bb_pval(allele_biased_rep)*bonferroni_correction_factor);

    ceres_names = {pstring1; pstring2};	
    
    for i = 1:N_bones    
        ceres = NaN(L(i),1);
        ceres(samples{i}.allele == allele_equal_rep ) = 1;
        ceres(samples{i}.allele == allele_biased_rep) = 2;
        plot_allele_on_transcriptome('EB', samples{i}.y, -samples{i}.x, ceres, ceres_names, [3; 3], ...
                                     ['br']', ['oo']', [1; 2], true, [2; 1], title_str{i}, true, []);
        paper_print(sprintf('%s/Figures/EB/BoneBias%sWithLegend', results_dir, bones{i}));
        plot_allele_on_transcriptome('EB', samples{i}.y, -samples{i}.x, ceres, ceres_names, [3; 3], ...
                                     ['br']', ['oo']', [1; 2], false, [2; 1], title_str{i}, true, []);
        paper_print(sprintf('%s/Figures/EB/BoneBias%sWithoutLegend', results_dir, bones{i}));
    end
    
    hsc_pheno    = find(strcmp(coarse_grain_pheno(:,1), 'HSC'));
    
    ceres_names = {'Parent HSC'; 'Childless HSC'; 'HSC-Derived Progeny'; 'Progeny-Only Clone'};    
    
    ceres = NaN(size(pooled.x));    
    
    ceres(pooled.cg_pheno == hsc_pheno & ismember(pooled.allele, combined.sig_alleles.hsc_derived_allele )) = 1;
    ceres(pooled.cg_pheno == hsc_pheno & ismember(pooled.allele, combined.sig_alleles.hsc_only_allele    )) = 2;
    ceres(pooled.cg_pheno ~= hsc_pheno & ismember(pooled.allele, combined.sig_alleles.hsc_derived_allele )) = 3;
    ceres(pooled.cg_pheno ~= hsc_pheno & ismember(pooled.allele, combined.sig_alleles.derived_only_allele)) = 4;
        
    color = [0 0.7 0;
             0 0 1;
             0 0.7 0;
             0.7 0.7 0];
          
    marker = ['oodo']';
          
    highlight_patch = get_highlight_patch(pooled.x, pooled.y, find(pooled.cg_pheno==hsc_pheno));
    highlight_patch = [pooled.y(highlight_patch) -pooled.x(highlight_patch)];
              
    temp_title = sprintf('All Bones - %d sig. clones with %d cells', ...
        length(combined.sig_alleles.all), sum(combined.summary.allele_freqs(combined.sig_alleles.all)));
    
    plot_allele_on_transcriptome('EB', pooled.y, -pooled.x, ceres, ceres_names, [3; 12; 3; 12], ...
                                 color, marker, [2; 4; 1; 3], true, [1; 3; 2; 4], temp_title, false, highlight_patch);
    paper_print(sprintf('%s/Figures/EB/ClonalExpansionPooledWithLegend', results_dir));
    
    plot_allele_on_transcriptome('EB', pooled.y, -pooled.x, ceres, ceres_names, [3; 12; 3; 12], ...
                                 color, marker, [2; 4; 1; 3], false, [1; 3; 2; 4], temp_title, false, highlight_patch);
    paper_print(sprintf('%s/Figures/EB/ClonalExpansionPooledWithoutLegend', results_dir));
      
    plot_EB_clone_bone_pheno_heatmap(combined.clone_by_cg_pheno_by_bone_breakdown(candidate_bb_alleles,:,:), ...
                                     candidate_bb_alleles, bb_alleles, bb_pval(bb_alleles)*bonferroni_correction_factor, ...
                                     combined.hsc_derived_allele, combined.hsc_only_allele, combined.derived_only_allele, color([1 2 4],:), ...
                                     coarse_grain_pheno(:,1), title_str);
    paper_print(sprintf('%s/Figures/EB/Heatmap', results_dir)); 
  
    close all;
    
end
   
function pstring = get_pstring(allele, pval)
    pval = min(pval, 1.0);
    if (pval < 1e-6)
        pstring = sprintf('Clone %d (p < 1e-6)', allele);
    elseif (pval >= 0.1)
        pstring = sprintf('Clone %d (p = %.1g)', allele, pval);
    else
        pstring = sprintf('Clone %d (p = %2.1e)', allele, pval);
    end
end