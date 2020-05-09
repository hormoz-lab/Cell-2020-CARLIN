function generate_SC_stats(results_dir)

    load(sprintf('%s/5FU/Analysis.mat', results_dir));
    
    FU5_ctrl = compute_conditioned_clone_and_cell_populations(samples, hsc_pheno);
    ctrl = compute_conditioned_clone_and_cell_populations(samples(neg_samples), hsc_pheno);
    FU5  = compute_conditioned_clone_and_cell_populations(samples(pos_samples), hsc_pheno);
  
    [FU5.z_derived_is_hsc_derived, FU5.p_derived_is_hsc_derived] = ...
        two_proportion_ztest( FU5.cells.derived_in_hsc_derived,  FU5.cells.in_derived_only, ...
                             ctrl.cells.derived_in_hsc_derived, ctrl.cells.in_derived_only);
                         
    [FU5.z_clone_is_hsc_rooted, FU5.p_clone_is_hsc_rooted] = ...
        two_proportion_ztest( FU5.clones.hsc_rooted,  FU5.clones.sig- FU5.clones.hsc_rooted, ...
                             ctrl.clones.hsc_rooted, ctrl.clones.sig-ctrl.clones.hsc_rooted);
                         
    [FU5.z_cell_in_hsc_rooted, FU5.p_cell_in_hsc_rooted] = ...
        two_proportion_ztest( FU5.cells.in_hsc_rooted,  FU5.cells.in_sig- FU5.cells.in_hsc_rooted, ...
                             ctrl.cells.in_hsc_rooted, ctrl.cells.in_sig-ctrl.cells.in_hsc_rooted);
                         
    [FU5.z_clone_is_hsc_derived, FU5.p_clone_is_hsc_derived] = ...
        two_proportion_ztest( FU5.clones.hsc_derived,  FU5.clones.sig- FU5.clones.hsc_derived, ...
                             ctrl.clones.hsc_derived, ctrl.clones.sig-ctrl.clones.hsc_derived);
                         
    [FU5.z_cell_in_hsc_derived, FU5.p_cell_in_hsc_derived] = ...
        two_proportion_ztest( FU5.cells.in_hsc_derived,  FU5.cells.in_sig- FU5.cells.in_hsc_derived, ...
                             ctrl.cells.in_hsc_derived, ctrl.cells.in_sig-ctrl.cells.in_hsc_derived);
 
    FU5_ctrl.cells.hsc_childless = sum(cellfun(@(x) sum(x.cg_pheno==hsc_pheno & ismember(x.allele, x.sig_alleles.hsc_only_allele   )), samples));
    FU5_ctrl.cells.hsc_parent    = sum(cellfun(@(x) sum(x.cg_pheno==hsc_pheno & ismember(x.allele, x.sig_alleles.hsc_derived_allele)), samples));
                         
    FU5_ctrl.cells.parent_HSC_in_parent_cluster       = sum(cellfun(@(x) sum(ismember(x.louvain, parent_HSC_clusters)    & ismember(x.allele, x.sig_alleles.hsc_derived_allele)), samples));
    FU5_ctrl.cells.parent_HSC_in_childless_cluster    = sum(cellfun(@(x) sum(ismember(x.louvain, childless_HSC_clusters) & ismember(x.allele, x.sig_alleles.hsc_derived_allele)), samples));
    FU5_ctrl.cells.childless_HSC_in_parent_cluster    = sum(cellfun(@(x) sum(ismember(x.louvain, parent_HSC_clusters)    & ismember(x.allele, x.sig_alleles.hsc_only_allele   )), samples));
    FU5_ctrl.cells.childless_HSC_in_childless_cluster = sum(cellfun(@(x) sum(ismember(x.louvain, childless_HSC_clusters) & ismember(x.allele, x.sig_alleles.hsc_only_allele   )), samples));
    
    [FU5_ctrl.z_parent_HSC_in_parent_cluster, FU5_ctrl.p_parent_HSC_in_parent_cluster] = ...
        two_proportion_ztest(FU5_ctrl.cells.parent_HSC_in_parent_cluster,    FU5_ctrl.cells.childless_HSC_in_parent_cluster, ...
                             FU5_ctrl.cells.parent_HSC_in_childless_cluster, FU5_ctrl.cells.childless_HSC_in_childless_cluster);
                         
    prof_label = splitlines(fileread(sprintf('%s/5FU/Seurat/DGE_HSC_Prolif_Labels.csv', results_dir)));
    prof_label = prof_label(1:end-1);
    [is, where] = cellfun(@(x) ismember(x, {'"Low"'; '"Mid"'; '"High"'; '"Other"'}), prof_label);
    assert(all(is));
    prof_pop = accumarray([pooled.louvain+1, where],1);
    FU5_ctrl.cells.high_prolif_in_parent_cluster    = sum(prof_pop(parent_HSC_clusters+1, 3));
    FU5_ctrl.cells.low_prolif_in_parent_cluster     = sum(sum(prof_pop(parent_HSC_clusters+1, [1:2])));
    FU5_ctrl.cells.high_prolif_in_childless_cluster = sum(prof_pop(childless_HSC_clusters+1, 3));
    FU5_ctrl.cells.low_prolif_in_childless_cluster  = sum(sum(prof_pop(childless_HSC_clusters+1, [1:2])));
    
    [FU5_ctrl.z_high_prolif_in_parent_cluster, FU5_ctrl.p_high_prolif_in_parent_cluster] = ...
        two_proportion_ztest(FU5_ctrl.cells.high_prolif_in_parent_cluster, FU5_ctrl.cells.low_prolif_in_parent_cluster, ...
                             FU5_ctrl.cells.high_prolif_in_childless_cluster, FU5_ctrl.cells.low_prolif_in_childless_cluster);
    
    FU5_ctrl.cells.high_prolif_in_parent_HSC    = sum(cellfun(@(x,y) sum(ismember(x.allele, x.sig_alleles.hsc_derived_allele) & strcmp(y,'"High"')), samples, mat2cell(prof_label, accumarray(pooled.sample,1))));
    FU5_ctrl.cells.high_prolif_in_childless_HSC = sum(cellfun(@(x,y) sum(ismember(x.allele, x.sig_alleles.hsc_only_allele) & strcmp(y,'"High"')), samples, mat2cell(prof_label, accumarray(pooled.sample,1))));                   
    
    [FU5_ctrl.z_high_prolif_in_parent_HSC, FU5_ctrl.p_high_prolif_in_parent_HSC] = ...
        two_proportion_ztest(FU5_ctrl.cells.high_prolif_in_parent_HSC, FU5_ctrl.cells.hsc_parent-FU5_ctrl.cells.high_prolif_in_parent_HSC, ...
                             FU5_ctrl.cells.high_prolif_in_childless_HSC, FU5_ctrl.cells.hsc_childless-FU5_ctrl.cells.high_prolif_in_childless_HSC);

    
    rng(43);
    for i = 1:length(pos_samples)
        for j = 1:length(neg_samples)
            FU5.nc{i,j} = compare_number_of_clones(samples{pos_samples(i)}, samples{neg_samples(j)});
        end
    end
                         
    rng(1034309243);
    for i = 1:length(pos_samples)
        for j = 1:length(neg_samples)
            FU5.acs{i,j} = compare_clone_size_distribution(samples{pos_samples(i)}, samples{neg_samples(j)});            
        end        
    end
    
    rng(5690934);
    for i = 1:length(pos_samples)
        FU5.csd{i} = compare_clone_size_distribution(samples{pos_samples(i)});
    end
    
    for i = 1:length(neg_samples)
        FU5.csd{length(pos_samples)+i} = compare_clone_size_distribution(samples{neg_samples(i)});
    end
    
    clearvars -except ctrl FU5 FU5_ctrl results_dir samples neg_samples

    load(sprintf('%s/EB/Analysis.mat', results_dir), 'combined', 'pooled', 'hsc_pheno');
    
    EB.cells.derived_in_hsc_derived = sum(pooled.cg_pheno~=hsc_pheno & ismember(pooled.allele, combined.sig_alleles.hsc_derived_allele));
    EB.cells.in_derived_only = sum(ismember(pooled.allele, combined.sig_alleles.derived_only_allele));
        
    [EB.z_derived_in_hsc_derived, EB.p_derived_in_hsc_derived] = ...
        two_proportion_ztest(  EB.cells.derived_in_hsc_derived,   EB.cells.in_derived_only, ...
                             ctrl.cells.derived_in_hsc_derived, ctrl.cells.in_derived_only);
                         
    rng(445343);    
    for i = 1:length(neg_samples)
        EB.nc{i} = compare_number_of_clones(combined, samples{neg_samples(i)});
    end
            
                         
    rng(1383896);
    for i = 1:length(neg_samples)
        EB.acs{i} = compare_clone_size_distribution(combined, samples{neg_samples(i)});
    end
    
    rng(84854);
    EB.csd = compare_clone_size_distribution(combined);
    
    outdir = sprintf('%s/SC', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    save(sprintf('%s/Stats.mat', outdir), 'FU5', 'EB', 'ctrl', 'FU5_ctrl');

end

function out = compute_conditioned_clone_and_cell_populations(samples, hsc_pheno)

    out.clones.sig          = sum(cellfun(@(x) length(x.sig_alleles.all), samples));
    out.clones.hsc_derived  = sum(cellfun(@(x) length(x.sig_alleles.hsc_derived_allele), samples));
    out.clones.hsc_only     = sum(cellfun(@(x) length(x.sig_alleles.hsc_only_allele), samples));    
    out.clones.derived_only = sum(cellfun(@(x) length(x.sig_alleles.derived_only_allele), samples));
    out.clones.hsc_rooted   = out.clones.hsc_derived + out.clones.hsc_only;
    
    out.cells.in_sig          = sum(cellfun(@(x) sum(ismember(x.allele, x.sig_alleles.all)), samples));
    out.cells.in_hsc_derived  = sum(cellfun(@(x) sum(ismember(x.allele, x.sig_alleles.hsc_derived_allele)), samples));
    out.cells.in_hsc_only     = sum(cellfun(@(x) sum(ismember(x.allele, x.sig_alleles.hsc_only_allele)), samples));
    out.cells.in_derived_only = sum(cellfun(@(x) sum(ismember(x.allele, x.sig_alleles.derived_only_allele)), samples));
    out.cells.in_hsc_rooted   = out.cells.in_hsc_derived + out.cells.in_hsc_only;
    
    out.cells.hscs_in_hsc_derived    = sum(cellfun(@(x) sum(x.cg_pheno == hsc_pheno & ismember(x.allele, x.sig_alleles.hsc_derived_allele)), samples));
    out.cells.derived_in_hsc_derived = sum(cellfun(@(x) sum(x.cg_pheno ~= hsc_pheno & ismember(x.allele, x.sig_alleles.hsc_derived_allele)), samples));
    
end
