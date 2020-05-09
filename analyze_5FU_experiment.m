function analyze_5FU_experiment(analysis_dir, results_dir)

    mouse = {'FO898'; 'SB133'; 'FO892_1'; 'FO892_2'; 'FO897'};    
    data_path = {sprintf('%s/5FU/FO898/SC', analysis_dir);
                 sprintf('%s/5FU/SB133/SC', analysis_dir);  
                 sprintf('%s/5FU/FO892/SC/1', analysis_dir);
                 sprintf('%s/5FU/FO892/SC/2', analysis_dir);
                 sprintf('%s/5FU/FO897/SC', analysis_dir);
                 };
             
    N_mice = length(mouse);
    samples = cell(N_mice,1);
    
    pos_samples = [2 1];
    neg_samples = [4 3 5];
    
    
    % Annotations from experimentalists based on Seurat dot plots
    % and feature plots of marker genes    
    phenotypes = {'HSC',    [0 2 17 19];
                  'MPP',    [1 3 4 5 11 13 15 20 21 23 24 26 28 30];
                  'MyP',    [10 27];
                  'Neu',    [8 16];
                  'Mono',   [7 22];            
                  'DC',     [32];
                  'Ba',     [25];           
                  'Eos',    [];             
                  'Mk',     [9 18 31];
                  'Ery',    [6 12 14];       
                  'Ly',     [33];
                  'Other',  [29]};
        
    parent_HSC_clusters = [0 2];
    childless_HSC_clusters = [17 19];
    
    assert(isequal(sort(horzcat(phenotypes{:,2})), [0:max(horzcat(phenotypes{:,2}))]))
              
    coarse_grain_pheno = {'HSC',   [1];
                          'MPP',   [2];
                          'My',    [3:8]; 
                          'Mk',    [9];
                          'Ery',   [10];
                          'Ly',    [11];
                          'Other', [12]};
                      
    hsc_pheno = strcmp(coarse_grain_pheno(:,1), 'HSC');   
                      
    [folder, ~, ~] = fileparts(mfilename('fullpath'));    
    marker_genes = splitlines(fileread(sprintf('%s/marker_genes.txt', folder)));   
    clear folder;
    
    dp_data = parse_seurat_dotplot_output(sprintf('%s/5FU/Seurat', results_dir), marker_genes, horzcat(phenotypes{:,2})');
 
    for i = 1:N_mice
        samples{i} = parse_seurat_cellid_output(data_path{i});
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
    
    fn = fieldnames(samples{1});
    pooled = cellfun(@(x) struct2cell(x)', samples, 'UniformOutput', false);
    pooled = vertcat(pooled{:});
    pooled = arrayfun(@(i) vertcat(pooled{:,i}), [1:length(fn)], 'un', false);
    pooled = cell2struct(pooled, fn, 2);
    pooled.sample = repelem([1:length(samples)]', cellfun(@(x) length(x.allele), samples));
    pooled = rmfield(pooled, 'allele');
    clear fn;
    
    load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir));
    fdr_level = 0.05;
    
    for i = 1:N_mice
        
        samples{i}.clone_by_fg_pheno_breakdown = accumarray([samples{i}.allele+1 samples{i}.fg_pheno], 1, [max(samples{i}.allele)+1, size(phenotypes,1)]);
        samples{i}.clone_by_fg_pheno_breakdown = samples{i}.clone_by_fg_pheno_breakdown(2:end,:);
        assert(all(any(samples{i}.clone_by_fg_pheno_breakdown,2)));
        
        samples{i}.clone_by_cg_pheno_breakdown = accumarray([samples{i}.allele+1 samples{i}.cg_pheno], 1, [max(samples{i}.allele)+1, size(coarse_grain_pheno,1)]);
        samples{i}.clone_by_cg_pheno_breakdown = samples{i}.clone_by_cg_pheno_breakdown(2:end,:);
        assert(all(any(samples{i}.clone_by_cg_pheno_breakdown,2)));        
        
        samples{i}.fg_fb_pval = fate_bias_pvalue(samples{i}.summary, samples{i}.clone_by_fg_pheno_breakdown);
        samples{i}.cg_fb_pval = fate_bias_pvalue(samples{i}.summary, samples{i}.clone_by_cg_pheno_breakdown);
        
        samples{i}.clonal_pval.all = bank.compute_clonal_pvalue(samples{i}.summary);
        p_crit = benjamini_hochberg(samples{i}.clonal_pval.all(2:end), fdr_level);
        samples{i}.sig_alleles.all = find(samples{i}.clonal_pval.all <= p_crit);
                
        hsc_mask     = logical(samples{i}.clone_by_cg_pheno_breakdown(:,hsc_pheno));
        derived_mask = any(samples{i}.clone_by_cg_pheno_breakdown(:,~hsc_pheno),2);
        
        samples{i}.hsc_only_allele     = find( hsc_mask & ~derived_mask);
        samples{i}.hsc_derived_allele  = find( hsc_mask &  derived_mask);
        samples{i}.derived_only_allele = find(~hsc_mask &  derived_mask);        
        
        samples{i}.sig_alleles.hsc_only_allele     = intersect(samples{i}.sig_alleles.all, samples{i}.hsc_only_allele    );
        samples{i}.sig_alleles.hsc_derived_allele  = intersect(samples{i}.sig_alleles.all, samples{i}.hsc_derived_allele );
        samples{i}.sig_alleles.derived_only_allele = intersect(samples{i}.sig_alleles.all, samples{i}.derived_only_allele);
        
    end
    
    clear i derived_mask hsc_mask bank;
    
    hsc_pheno = find(hsc_pheno);
    
    create_seurat_HSC_type_labels(samples, hsc_pheno, results_dir);
    
    create_seurat_HSC_cluster_labels(samples, parent_HSC_clusters, childless_HSC_clusters, results_dir);

    outdir = sprintf('%s/5FU', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    clearvars analysis_dir results_dir
    save(sprintf('%s/Analysis.mat', outdir));
   
end

function out = parse_seurat_cellid_output(analysis_dir)

    load(sprintf('%s/Amplicon/Summary.mat', analysis_dir));
    umap = readtable(sprintf('%s/Transcriptome/Seurat/UMapCoordinates.csv', analysis_dir));    
    assert(isequal(ref_CBs, cellfun(@(x) x(1:find(x=='-')-1), umap.Var1, 'un', false)));
    
    out.x = umap.Var2;
    out.y = umap.Var3;
    out.summary = summary;
    
    louvain = readtable(sprintf('%s/Transcriptome/Seurat/LouvainClusters.csv', analysis_dir));
    out.louvain = cellfun(@str2double, louvain.Var2);

    out.allele = NaN(size(out.louvain));        
    
    for j = 1:size(summary.alleles,1)
        [is, which_CB] = ismember(summary.allele_colony{j}, ref_CBs);
        assert(all(is));
        out.allele(which_CB) = j;
    end
    
end

function create_seurat_HSC_type_labels(samples, hsc_pheno, results_dir)

    label = cellfun(@(x) repmat({'Other'}, [length(x.louvain), 1]), samples, 'un', false);
    for i = 1:length(samples)
        label{i}(samples{i}.cg_pheno == hsc_pheno & ismember(samples{i}.allele, samples{i}.sig_alleles.hsc_only_allele)) = {'Childless'};
        label{i}(samples{i}.cg_pheno == hsc_pheno & ismember(samples{i}.allele, samples{i}.sig_alleles.hsc_derived_allele)) = {'Parent'};
    end
        
    label = vertcat(label{:});
    
    outdir = sprintf('%s/5FU/Seurat/', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    fid = fopen(sprintf('%s/DGE_HSC_Type_Labels.csv', outdir), 'wt');
    fprintf(fid, '%s\n', label{:});
    fclose(fid);
    
end


function create_seurat_HSC_cluster_labels(samples, parent_clusters, childless_clusters, results_dir)

    label = cellfun(@(x) repmat({'Other'}, [length(x.louvain), 1]), samples, 'un', false);
    for i = 1:length(samples)
        label{i}(ismember(samples{i}.louvain, childless_clusters)) = {'Childless'};
        label{i}(ismember(samples{i}.louvain, parent_clusters))    = {'Parent'};
    end
    
    label = vertcat(label{:});
    
    outdir = sprintf('%s/5FU/Seurat', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    fid = fopen(sprintf('%s/DGE_HSC_Cluster_Labels.csv', outdir), 'wt');
    fprintf(fid, '%s\n', label{:});
    fclose(fid);

end
