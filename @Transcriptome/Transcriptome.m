% A class to help with preprocessing and storage of count matrices. There's
% nothing new here. I found it convenient to have something in MATLAB where
% I could also play around with things, but you can use your favorite
% package instead for working with transcriptomes. CARLIN just requires a
% reference barcode list for denoising amplicon barcodes. There are also
% some helper functions for preparing inputs to the SPRING GUI, but this is
% pretty cursory. All this effectively does for the analysis is pick a UMI 
% cutoff and produces a reference barcode list where cells have at least that
% many UMIs and low mitochondrial fraction.

classdef Transcriptome
    
    properties (SetAccess = immutable, GetAccess = public)
        source_file
        barcodes
        genes
        counts
        normalized
        scaled
        cutoffs
        masks
        platform
        metrics
        params
    end
    
        
    methods (Static)
        
        % Loading SC count matrices
        [counts, barcodes, genes, platform] = load_count_matrix(input_count_file);
        [counts, barcodes, genes] = load_indrops_matrix(count_file);
        [counts, barcodes, features] = load_10x_h5_matrix(h5_file);
        
        % Normalizations
        scaled = normalize_genes(normalized, recentre);
        normalized = normalize_cells(counts, exclude_frac);
        
        % Filters
        missing_gene = filter_out_genes_missing(A);
        mt_mask = filter_out_genes_mt(genes);
        [var_masks, vscores] = filter_out_genes_low_variability(counts, min_mean, min_count, min_cell, v_pctl);
        [cc_mask, cc_uncorr_mask, rho, pval, cc_signature] = filter_out_genes_cc_corr(normalized, scaled, genes, corr_cutoff)

        [mt_mask, mt_frac_cutoff] = filter_out_cells_high_mt_frac(A, mt_gene_mask, mt_frac_cutoff);        
        [umi_mask, umi_cutoff] = filter_out_cells_low_umi(A, umi_cutoff);
        
         % Some helper utilities for InDrops matrices
        remap_indrops_internal_bc(int_bc_file, abundant_bc_file, outfile);
        
        % Perform standard QC checks
        prepare_QC_plot(counts, umi_mt_cell_mask, mt_genes, cutoffs);
            
        % Merging transcriptomes - really only using to get rid of the same
        % cell cycle genes in 
        function [summary, sample_map, allele_breakdown_by_sample] = prepare_SPRING_input(transcriptome, amplicon, outdir, summary_names)
            
            Nfiles = size(transcriptome,1);
            if (Nfiles > 1)
                assert(size(amplicon,1) == Nfiles);
                assert(nargin == 4);
                assert(size(summary_names,1) == Nfiles);
            else
                transcriptome = {transcriptome};                
                amplicon = {amplicon};
            end
            
            gene_mask = [1:length(transcriptome{1}.genes)]';
            for i = 1:Nfiles
                gene_mask = intersect(gene_mask, transcriptome{1}.masks.genes.cc_uncorr_mask);
            end
            
            if (~exist(outdir, 'dir'))
                mkdir(outdir);        
            end
            
            outfile = sprintf('%s/filtered_genes_cc.txt', outdir);
            fprintf('Writing gene list file: %s\n', outfile);
            fid = fopen(outfile, 'wt');
            fprintf(fid, '%s\n', transcriptome{1}.genes{gene_mask(1:end-1)});
            fprintf(fid, '%s', transcriptome{1}.genes{gene_mask(end)});
            fclose(fid);
            
            summary = vertcat(amplicon{:});
            [summary, sample_map, allele_breakdown_by_sample] = ExperimentSummary.FromMerge(vertcat(summary.summary));
            mark_singleton = size(allele_breakdown_by_sample,1) > 499;

            cluster_labels = cell(Nfiles,1);

            for i = 1:Nfiles
                assert(isequal(amplicon{i}.ref_CBs, transcriptome{i}.barcodes(transcriptome{i}.masks.cells.joint)));
                cluster_labels{i} = repmat(cellstr('Allele?'), [size(amplicon{i}.ref_CBs,1),1]);
                for j = 1:size(amplicon{i}.summary.alleles,1)
                    [is, which_CB] = ismember(amplicon{i}.summary.allele_colony{j}, amplicon{i}.ref_CBs);
                    assert(all(is));
                    if (mark_singleton && summary.allele_freqs(sample_map{i}(j))==1)
                        cluster_labels{i}(which_CB) = {sprintf('AlleleSingleton')};
                    else
                        cluster_labels{i}(which_CB) = {sprintf('Allele%d', sample_map{i}(j))};
                    end
                end    
                assert(sum(~strcmp(cluster_labels{i}, 'Allele?')) == amplicon{i}.summary.N.called_tags);
            end
            
            counts = cellfun(@(x) x.counts(x.masks.cells.joint, gene_mask), transcriptome, 'un', false);
            
            joint_cluster_labels = vertcat(cluster_labels{:});
            
            outfile = sprintf('%s/CellLabels.txt', outdir);
            fprintf('Writing cell labels file: %s\n', outfile);
            fid = fopen(outfile, 'wt');
            if (Nfiles > 1)
                fprintf(fid, 'Bone');
                bone_labels = repelem(summary_names, cellfun(@(x) size(x,1), counts));
                fprintf(fid, ',%s', bone_labels{:});
                fprintf(fid, '\n');
            end
            fprintf(fid, 'CARLIN Allele');
            fprintf(fid, ',%s', joint_cluster_labels{:});
            fclose(fid);
            
            outfile = sprintf('%s/filtered_cells_umi_mt_genes_cc.tsv', outdir);
            fprintf('Writing TSV file: %s\n', outfile);
            dlmwrite(outfile, full(vertcat(counts{:})), 'Delimiter', '\t');
            fprintf('Zipping TSV file: %s\n', outfile);
            gzip(outfile);
            delete(outfile);           
        end
        
        function spring_data = get_spring_output(spring_path)
            
            assert(exist(spring_path, 'dir')==7);
            coords = cellfun(@(x) strsplit(x, ','), splitlines(fileread(sprintf('%s/coordinates.csv', spring_path))), 'un', false);
            coords = str2double(vertcat(coords{1:end-1}));
            spring_data.coords = coords(:,2:3);

            labels = cellfun(@(x) strsplit(x, ','), splitlines(fileread(sprintf('%s/cell_groupings.csv', spring_path))), 'un', false);
            labels = vertcat(labels{1:end-1})';

            spring_data.louvain = str2double(labels(2:end, strcmp(labels(1,:), 'Louvain cluster')));

            allele = labels(2:end, strcmp(labels(1,:), 'CARLIN Allele'));
            allele = strrep(allele, 'Allele', '');
            allele = strrep(allele, '?', 'NaN');
            allele = strrep(allele, 'Singleton', 'inf');
            spring_data.allele = str2double(allele);

            if (any(strcmp(labels(1,:), 'Bone')))
                spring_data.sample_names = labels(2:end, strcmp(labels(1,:), 'Bone'));
            end
        end
    end
    
    methods (Access = public)
        
        function obj = Transcriptome(input_count_file, params, outdir)
            
            if (nargin < 2)
                params = struct;
            end
            
            obj.source_file = input_count_file;
            [obj.counts, obj.barcodes, obj.genes, obj.platform] = Transcriptome.load_count_matrix(input_count_file);
            
            % Exclude dominant genes for cell normalization
            if (~isfield(params, 'exclude_frac'))
                params.exclude_frac = 1.0;
            end            
            obj.normalized = Transcriptome.normalize_cells(obj.counts, params.exclude_frac);
                        
            if (~isfield(params, 'normalization_recentre'))
                params.normalization_recentre = true;
            end            
            obj.scaled = Transcriptome.normalize_genes(obj.normalized, params.normalization_recentre);
            
            if (nargin == 3)
                if (~exist(outdir, 'dir'))
                    mkdir(outdir);        
                end
            end
            
            if (~isfield(params, 'min_mean'))
                params.min_mean = 0;
            end
            if (~isfield(params, 'min_count'))
                params.min_count = 3;
            end
            if (~isfield(params, 'min_cell'))
                params.min_cell = 5;
            end
            if (~isfield(params, 'v_pctl'))
                params.v_pctl = 85;
            end
            
            [obj.masks.genes, obj.metrics.genes.vscores] = Transcriptome.filter_out_genes_low_variability(...
                obj.counts, params.min_mean, params.min_count, params.min_cell, params.v_pctl);
            
            if (nargin == 3)
                print(sprintf('%s/GeneVariability.png', outdir),'-dpng','-r0');
                close;
            end
            
            obj.masks.genes.missing = Transcriptome.filter_out_genes_missing(obj.counts);
            obj.masks.genes.mt      = Transcriptome.filter_out_genes_mt(obj.genes);
            
            if (~isfield(params, 'umi_cutoff'))
                params.umi_cutoff = -1;
            end            
            [obj.masks.cells.umi, obj.cutoffs.umi] = Transcriptome.filter_out_cells_low_umi(obj.counts, params.umi_cutoff);
            
            if (~isfield(params, 'mt_frac_cutoff'))
                params.mt_frac_cutoff = 0.15;
            end
            [obj.masks.cells.mt, obj.cutoffs.mt_frac] = Transcriptome.filter_out_cells_high_mt_frac(...
                obj.counts, obj.masks.genes.mt, params.mt_frac_cutoff);
            
            obj.masks.cells.joint = obj.masks.cells.mt & obj.masks.cells.umi;
            fprintf('Cells surviving joint filter: %d/%d\n', sum(obj.masks.cells.joint), length(obj.masks.cells.joint));
            
            Transcriptome.prepare_QC_plot(obj.counts, obj.masks.cells, obj.masks.genes.mt, obj.cutoffs);
            obj.masks.cells.mt    = uint16(find(obj.masks.cells.mt));
            obj.masks.cells.umi   = uint16(find(obj.masks.cells.umi));
            obj.masks.cells.joint = uint16(find(obj.masks.cells.joint));
            
            if (nargin == 3)
                print(sprintf('%s/CellQC.png', outdir),'-dpng','-r0');
                close;
            end
            
            if (~isfield(params, 'corr_cutoff'))
                params.corr_cutoff = 0.1;
            end
            
            [obj.masks.genes.cc_mask, obj.masks.genes.cc_uncorr_mask, obj.metrics.genes.cc_corr_rho, ...
             obj.metrics.genes.cc_corr_pval, obj.metrics.genes.cc_signature] = ...
                Transcriptome.filter_out_genes_cc_corr(obj.normalized, obj.scaled, obj.genes, params.corr_cutoff);
            
            obj.params = params;
                        
        end
        
        function prepare_CARLIN_input(obj, outdir)
        
            if (~exist(outdir, 'dir'))
                mkdir(outdir);        
            end
            if (strcmp(obj.platform, '10X'))
                outfile = sprintf('%s/filtered_barcodes_umi_mt.txt', outdir);
            else
                outfile = sprintf('%s/filtered_intbcs_umi_mt.txt', outdir);
            end
            fprintf('Writing barcode file: %s\n', outfile);
            fid = fopen(outfile, 'wt');
            fprintf(fid, '%s\n', obj.barcodes{obj.masks.cells.joint(1:end-1)});
            fprintf(fid, '%s', obj.barcodes{obj.masks.cells.joint(end)});
            fclose(fid);
            
        end
    end
end
