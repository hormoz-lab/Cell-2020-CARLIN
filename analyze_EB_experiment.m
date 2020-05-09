function analyze_EB_experiment(analysis_dir, results_dir)

    load(sprintf('%s/EB/FO864/SC/Combined/Amplicon/Summary.mat', analysis_dir));    
    
    bones = strrep(sample_names, 'FO864/', '');
    N_bones = length(bones);
    
    % Annotations from experimentalists based on Seurat dot plots
    % and feature plots of marker genes
    phenotypes = {'HSC',    [0 7 13 28];
                  'MPP',    [1 2 8 9 12 19 21 22 25 30];
                  'MyP',    [15 16];
                  'Neu',    [3 11 29];
                  'Mono',   [5 6 17 35];
                  'DC',     [26];
                  'Ba',     [23];
                  'Eos',    [32];
                  'Mk',     [4 20 24 27];
                  'Ery',    [10 14 18 34];
                  'Ly',     [33];
                  'Other'   [31]};
              
    assert(isequal(sort(horzcat(phenotypes{:,2})), [0:max(horzcat(phenotypes{:,2}))]))
    
    coarse_grain_pheno = {'HSC',   [1];
                          'MPP',   [2];
                          'My',    [3:8]; 
                          'Mk',    [9];
                          'Ery',   [10];
                          'Ly',    [11];
                          'Other', [12]};
                      
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    marker_genes = splitlines(fileread(sprintf('%s/marker_genes.txt', folder)));
    clear folder;
    
    dp_data = parse_seurat_dotplot_output(sprintf('%s/EB/Seurat', results_dir), marker_genes, horzcat(phenotypes{:,2})');

    fdr_level = 0.05;

    for i = 1:N_bones
        samples{i} = parse_seurat_cellid_output(analysis_dir, combined, samples{i}, bones{i}, find(endsWith(sample_names, bones{i})));
        samples{i}.allele(isnan(samples{i}.allele)) = 0;
        
        samples{i}.fg_pheno = cellfun(@(x) ismember(samples{i}.louvain', x), phenotypes(:,2), 'un', false);
        samples{i}.fg_pheno = vertcat(samples{i}.fg_pheno{:});
        assert(all(sum(samples{i}.fg_pheno,1)==1));
        [samples{i}.fg_pheno, ~] = find(samples{i}.fg_pheno);
        
        samples{i}.cg_pheno = cellfun(@(x) ismember(samples{i}.fg_pheno', x), coarse_grain_pheno(:,2), 'un', false);
        samples{i}.cg_pheno = vertcat(samples{i}.cg_pheno{:});
        assert(all(sum(samples{i}.cg_pheno,1)==1));
        [samples{i}.cg_pheno, ~] = find(samples{i}.cg_pheno);   
    end
    
    clear i;
    
    fn = fieldnames(samples{1});
    pooled = cellfun(@(x) struct2cell(x)', samples, 'UniformOutput', false);
    pooled = vertcat(pooled{:});
    pooled = arrayfun(@(i) vertcat(pooled{:,i}), [1:length(fn)], 'un', false);
    pooled = cell2struct(pooled, fn, 2);
    L = cellfun(@(s) length(s.x), samples);
    pooled.sample = repelem([1:N_bones]', L);
    clear fn; 
    
    load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir));
    combined.clonal_pval.all = bank.compute_clonal_pvalue(combined.summary);
    clear bank;
    
    p_crit = benjamini_hochberg(combined.clonal_pval.all(2:end), fdr_level);
    combined.sig_alleles.all = find(combined.clonal_pval.all <= p_crit);
    
    sig_level = 0.05;
    
    [bb_pval, min_cells] = fate_bias_pvalue(combined.summary, combined.allele_breakdown_by_sample, sig_level);
    bones_present = sum(logical(combined.allele_breakdown_by_sample),2);

    nontrivial_pop_alleles       = find(combined.summary.allele_freqs >= min_cells);
    candidate_bb_alleles         = intersect(combined.sig_alleles.all, nontrivial_pop_alleles);    
    bonferroni_correction_factor = length(combined.sig_alleles.all);
    bb_alleles                   = intersect(candidate_bb_alleles, find(bb_pval <= sig_level/bonferroni_correction_factor));

    allele_biased_rep = intersect(bb_alleles, find(bones_present < N_bones));
    allele_biased_rep = allele_biased_rep(1);
    allele_equal_rep  = setdiff(candidate_bb_alleles, bb_alleles);
    allele_equal_rep  = allele_equal_rep(find(allele_equal_rep > allele_biased_rep, 1, 'first'));
            
    combined.clone_by_cg_pheno_by_bone_breakdown = accumarray([pooled.allele+1 pooled.cg_pheno pooled.sample], 1, ...
                                                     [max(pooled.allele)+1, length(coarse_grain_pheno), N_bones]);
    combined.clone_by_cg_pheno_by_bone_breakdown = combined.clone_by_cg_pheno_by_bone_breakdown(2:end,:,:);    
    assert(isequal(squeeze(sum(combined.clone_by_cg_pheno_by_bone_breakdown,2)), combined.allele_breakdown_by_sample));
    
    clone_by_cg_pheno = squeeze(sum(combined.clone_by_cg_pheno_by_bone_breakdown, 3));
    assert(all(any(clone_by_cg_pheno,2)));
    
    hsc_pheno = strcmp(coarse_grain_pheno(:,1), 'HSC');
        
    hsc_mask     = logical(clone_by_cg_pheno(:,hsc_pheno));
    derived_mask = any(clone_by_cg_pheno(:,~hsc_pheno),2);

    combined.hsc_only_allele     = find( hsc_mask & ~derived_mask);
    combined.hsc_derived_allele  = find( hsc_mask &  derived_mask);
    combined.derived_only_allele = find(~hsc_mask &  derived_mask);
    
    combined.sig_alleles.hsc_only_allele     = intersect(combined.sig_alleles.all, combined.hsc_only_allele    );
    combined.sig_alleles.hsc_derived_allele  = intersect(combined.sig_alleles.all, combined.hsc_derived_allele );
    combined.sig_alleles.derived_only_allele = intersect(combined.sig_alleles.all, combined.derived_only_allele);
    
    hsc_pheno = find(hsc_pheno);
    
    outdir = sprintf('%s/EB', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    clearvars analysis_dir results_dir
    save(sprintf('%s/Analysis.mat', outdir));
    
end
    
function out = parse_seurat_cellid_output(analysis_dir, combined, sample, which_bone, which_sample)

    umap = readtable(sprintf('%s/EB/FO864/SC/%s/Transcriptome/Seurat/UMapCoordinates.csv', analysis_dir, which_bone));    
    assert(isequal(sample.ref_CBs, cellfun(@(x) x(1:find(x=='-')-1), umap.Var1, 'un', false)));
    
    out.x = umap.Var2;
    out.y = umap.Var3;
    out.summary = sample.summary;
    
    % There's an outlier island that introduces a lot of interstitial white space
    out.x(out.x>12) = out.x(out.x>12)-6;
    
    louvain = readtable(sprintf('%s/EB/FO864/SC/%s/Transcriptome/Seurat/LouvainClusters.csv', analysis_dir, which_bone));
    out.louvain = cellfun(@str2double, louvain.Var2);
    
    out.allele = NaN(size(out.louvain));    
    
    for j = 1:size(out.summary.alleles,1)
        [is, which_CB] = ismember(out.summary.allele_colony{j}, sample.ref_CBs);
        assert(all(is));
        out.allele(which_CB) = combined.sample_map{which_sample}(j);
    end
end