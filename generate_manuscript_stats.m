function generate_manuscript_stats(processed_path, analysis_path, results_path)
    
    fid = fopen(sprintf('%s/ManuscriptStats.txt', results_path), 'wt');
    
    %% ES Cell Stats
    
    fprintf(fid, 'ES CELLS\n');
    
    load(sprintf('%s/ChronicInduction/PosDox/Low/96h/Summary.mat', analysis_path), 'summary');
    
    mut_events = cellfun(@(x) Mutation.identify_Cas9_events(x), summary.alleles, 'un', false);
    mut_events = vertcat(mut_events{:});
    del_events = mut_events(ismember(vertcat(mut_events.type), 'CD'));
    ins_events = mut_events(ismember(vertcat(mut_events.type), 'CI'));
    L_del = cellfun(@(x) length(degap(x)), {del_events.seq_old}');
    L_ins = cellfun(@(x) length(degap(x)), {ins_events.seq_new}');
    fprintf(fid, '\nMin deletion length: %d\n', min(L_del));
    fprintf(fid, '\nMax deletion length: %d\n', max(L_del));
    fprintf(fid, '\nMax insertion length: %d\n', max(L_ins));
    
    fprintf(fid, '\nLow Dox @ 96h: %d edited allles in %d edited cells with %d singletons\n', ...
        length(summary.alleles)-1, summary.N.eventful_tags, sum(summary.allele_freqs==1));

    %% In Vitro Phylogeny
    
    fprintf(fid, '\nIN VITRO PHYLOGENY\n');
    
    load(sprintf('%s/Trees/InVitroPhylogeny/Analysis.mat', results_path), 'fn', 'fp')    
    fprintf(fid, '\nAverage FP/FN rates (%%) among transcripts: [%8.6g %8.6g]\n', mean(fp)*100, mean(fn)*100);
    
    load(sprintf('%s/Trees/InVitroPhylogeny/Lineage.mat', results_path), 'allele_freq', 'N_stable_alleles', 'stable_allele_at_node');
    consensus = sum(allele_freq(stable_allele_at_node(2:N_stable_alleles)));
    fprintf(fid, '\nConsensus tree accounts for: %d of %d cells (%8.6g%%)\n', ...
        consensus, sum(allele_freq), consensus/sum(allele_freq)*100);

     
    %% Tissue Panel
     
    fprintf(fid, '\nTISSUE PANEL\n');
     
    [sample_sheet, samples_to_run] = get_experiment_metadata('TissuePanel');
    
    samples_to_run_pos = samples_to_run(contains(sample_sheet.ProcessedPath(samples_to_run), 'FO912'));
    samples_to_run_pos = samples_to_run_pos(~contains(sample_sheet.ProcessedPath(samples_to_run_pos), {'Brain'; 'Heart'; 'Muscle'}));
    summary = cell(length(samples_to_run_pos),1);
    for i = 1:length(samples_to_run_pos)
        [~, ~, output_path, ~] = prep_bulk_inputs(sample_sheet, samples_to_run_pos(i), processed_path, analysis_path);
        summary{i} = load(sprintf('%s/Summary.mat', output_path), 'summary');
    end
    summary = cellfun(@(x) x.summary, summary, 'un', false);
    pct_edited = cellfun(@(x) x.N.eventful_tags/x.N.called_tags, summary)*100;
    fprintf(fid, '\nRange in %% transcripts edited in +Dox tissue panel = [%8.6g, %8.6g]\n', min(pct_edited), max(pct_edited));
    
    samples_to_run_neg = samples_to_run(contains(sample_sheet.ProcessedPath(samples_to_run), 'FO916'));
    summary = cell(length(samples_to_run_neg),1);
    for i = 1:length(samples_to_run_neg)
        [~, ~, output_path, ~] = prep_bulk_inputs(sample_sheet, samples_to_run_neg(i), processed_path, analysis_path);
        summary{i} = load(sprintf('%s/Summary.mat', output_path), 'summary');
    end
    summary = cellfun(@(x) x.summary, summary, 'un', false);
    pct_edited = cellfun(@(x) x.N.eventful_tags/x.N.called_tags, summary)*100;
    fprintf(fid, '\nMean %% transcripts edited in -Dox tissue panel = %8.6g\n', mean(pct_edited));
    
    samples_to_run_pos = samples_to_run(contains(sample_sheet.ProcessedPath(samples_to_run), {'SB213'; 'SB214'; 'SB225'; 'SB226'}));    
    
    summary = cell(length(samples_to_run_pos),1);
    for i = 1:length(samples_to_run_pos)
        [~, ~, output_path, ~] = prep_bulk_inputs(sample_sheet, samples_to_run_pos(i), processed_path, analysis_path);
        summary{i} = load(sprintf('%s/Summary.mat', output_path), 'summary');
    end
    summary = cellfun(@(x) x.summary, summary, 'un', false);
    pct_edited = cellfun(@(x) x.N.eventful_tags/x.N.called_tags, summary)*100;
    fprintf(fid, '\nRange in %% transcripts edited in +Dox cell-sorted panel = [%8.6g, %8.6g]\n', min(pct_edited), max(pct_edited));
    
    %% Diversity
    
    fprintf(fid, '\nDIVERSITY\n');
    
    [sample_sheet, samples_to_run] = get_experiment_metadata('Diversity');
    samples_to_run = samples_to_run(contains(sample_sheet.ProcessedPath(samples_to_run), 'FO912'));
    
     
    pos_bank = load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_path));    
    is_template = ismember(cellfun(@(x) x.get_seq, pos_bank.bank.summary.alleles, 'un', false), ...
                           CARLIN_def.getInstance.seq.CARLIN);    
    pct_edited_by_mouse = sum(pos_bank.bank.allele_breakdown_by_sample(~is_template,:),1)./sum(pos_bank.bank.allele_breakdown_by_sample,1)*100;    
    fprintf(fid, '\nRange of %% edited transcripts per mouse: [%8.6g, %8.6g]\n', min(pct_edited_by_mouse), max(pct_edited_by_mouse));
    
    neg_bank = load(sprintf('%s/Banks/Protocol2/RNA/NegDox/Bank.mat', results_path));    
    is_template = ismember(cellfun(@(x) x.get_seq, neg_bank.bank.summary.alleles, 'un', false), ...
                           CARLIN_def.getInstance.seq.CARLIN);
    pct_edited = sum(neg_bank.bank.summary.allele_freqs(~is_template))/sum(neg_bank.bank.summary.allele_freqs)*100;
    fprintf(fid, '\n%% edited transcripts in -Dox mice: %8.6g\n', pct_edited);
    
    neg_cas9_bank = load(sprintf('%s/Banks/Protocol1/RNA/NegCas9/Bank.mat', results_path));
    is_template = ismember(cellfun(@(x) x.get_seq, neg_cas9_bank.bank.summary.alleles, 'un', false), ...
                           CARLIN_def.getInstance.seq.CARLIN);
    pct_edited = sum(neg_cas9_bank.bank.summary.allele_freqs(~is_template))/sum(neg_cas9_bank.bank.summary.allele_freqs)*100;
    fprintf(fid, '\n%% edited transcripts -Cas9 mice: %8.6g\n', pct_edited);

    unique_allele_mask = sum(logical(pos_bank.bank.allele_breakdown_by_sample),2)==1;
    unique_alleles_by_mouse = sum(logical(pos_bank.bank.allele_breakdown_by_sample(unique_allele_mask,:)),1);
    fprintf(fid, '\nRange Unique Alleles Per Mouse: [%d, %d]\n', min(unique_alleles_by_mouse), max(unique_alleles_by_mouse));
    pct_unique_alleles_by_mouse = unique_alleles_by_mouse./sum(logical(pos_bank.bank.allele_breakdown_by_sample),1)*100;
    fprintf(fid, '\nAverage %% Unique Alleles Per Mouse: %6.5g\n', mean(pct_unique_alleles_by_mouse));
    
    fprintf(fid, '\nDistinct Edited Alleles: %d\n', length(pos_bank.bank.summary.alleles)-1);
    
    fprintf(fid, '\nDistinct Edited Transcripts: %d\n', pos_bank.bank.model.transcripts.edited);
    
    fprintf(fid, '\nDiversity Estimate of Alleles: %d +/- %d\n', pos_bank.bank.model.alleles.estimated, round(pos_bank.bank.model.alleles.SE));    
    
    %% In Vivo Phylogeny
    
    fprintf(fid, '\nIN VIVO PHYLOGENY\n');
    
    load(sprintf('%s/Trees/InVivoPhylogeny/Lineage.mat', results_path), ...
        'allele_freq', 'N_stable_alleles', 'stable_allele_at_node', 'raw_combined');
    consensus = sum(allele_freq(stable_allele_at_node(2:N_stable_alleles)));
    total = sum(raw_combined.summary.allele_freqs(2:end));
    fprintf(fid, '\nConsensus tree accounts for: %d of %d edited transcripts (%8.6g%%)\n', ...
        consensus, total, consensus/total*100 );
    
    %% SC QC
    
    fprintf(fid, '\nSC QC\n');
    
    [sample_sheet, samples_to_run] = get_experiment_metadata({'5FU'; 'EB'});
    samples_to_run = samples_to_run(strcmp(sample_sheet.Type(samples_to_run), 'sc10xV3'));    
    summary = cell(length(samples_to_run),1);
    for i = 1:length(samples_to_run)
        [~, ~, output_path, ~] = prep_SC_inputs(sample_sheet, samples_to_run(i), processed_path, analysis_path);
        summary{i} = load(sprintf('%s/Summary.mat', output_path), 'summary');
    end
    summary = cellfun(@(x) x.summary, summary, 'un', false);
    CARLIN_capture = cellfun(@(x) x.N.called_tags/x.N.reference_tags*100, summary);    
    fprintf(fid, '\n[Min/Max CARLIN capture %%age: [%6.5g, %6.5g]\n', min(CARLIN_capture), max(CARLIN_capture));
    pct_edited = cellfun(@(x) x.N.eventful_tags/x.N.called_tags*100, summary);
    fprintf(fid, '\n[Min/Max edited CB %%age: [%6.5g, %6.5g]\n', min(pct_edited), max(pct_edited));
   
    [sample_sheet, samples_to_run] = get_experiment_metadata('SCReplicate');
    summary = cell(length(samples_to_run),1);
    for i = 1:length(samples_to_run)
        [~, ~, output_path, ~] = prep_SC_inputs(sample_sheet, samples_to_run(i), processed_path, analysis_path);
        summary{i} = load(sprintf('%s/Summary.mat', output_path), 'summary');
    end
    summary = cellfun(@(x) x.summary, summary, 'un', false);
    summary = reshape(summary, [2, 4])';
    cb_overlap = zeros(4,1);
    for i = 1:4
        [~, sample_map, allele_breakdown_by_sample] = ExperimentSummary.FromMerge([summary{i,1}; summary{i,2}]);
        allele_colony = cell(size(allele_breakdown_by_sample));
        for j = 1:2
            for k = 1:length(sample_map{j})
                allele_colony{sample_map{j}(k),j} = summary{i,j}.allele_colony{k};
            end
        end
        cb_overlap(i) = sum(arrayfun(@(k) length(intersect(allele_colony{k,1}, allele_colony{k,2})), [1:length(allele_colony)]')) / ...
                        sum(arrayfun(@(k) length(    union(allele_colony{k,1}, allele_colony{k,2})), [1:length(allele_colony)]'));
    end
    fprintf(fid, '\nSC Resampling overlap CBs: %6.5g +/- %6.5g\n', mean(cb_overlap), std(cb_overlap));
    
    
    %% EB
     
    fprintf(fid, '\nEB\n');

    load(sprintf('%s/EB/Analysis.mat', results_path));
    N_cells = cellfun(@(x) x.summary.N.reference_tags, samples);
    
    fprintf(fid, '\n[Min, max, total] number of cells in FO864 EB experiment = [%d, %d, %d]\n', min(N_cells), max(N_cells), sum(N_cells));
    
    fprintf(fid, '\nTotal number of Louvain clusters: %d\n', max(pooled.louvain)+1);
        
    fprintf(fid, '\nAll bones - Number of CARLIN alleles at FDR=%4.3g: %d/%d\n', ...
        fdr_level, length(combined.sig_alleles.all), length(combined.summary.alleles));
    num_sig_clones_per_bone = sum(logical(combined.allele_breakdown_by_sample(combined.sig_alleles.all,:)),1);    
    fprintf(fid, '\nAll bones - range of # significant clones across bones: [%d %d]\n', ...
        min(num_sig_clones_per_bone), max(num_sig_clones_per_bone));
    
    fprintf(fid, '\nAll bones - range of significant clone sizes: [%d %d]\n', ...
        min(combined.summary.allele_freqs(combined.sig_alleles.all)), ...
        max(combined.summary.allele_freqs(combined.sig_alleles.all)));
    
    progeny = length(combined.sig_alleles.hsc_derived_allele);
    total = length(union(combined.sig_alleles.hsc_derived_allele, combined.sig_alleles.hsc_only_allele));
    fprintf(fid, '\nAll bones - %d/%d (%8.6g%%) significant HSC-rooted clones have progeny\n', ...
        progeny, total, progeny/total*100);

    load(sprintf('%s/SC/Stats.mat', results_path), 'EB');
    
    non_HSC_in_clone = EB.cells.derived_in_hsc_derived;
    non_HSC = EB.cells.in_derived_only+EB.cells.derived_in_hsc_derived;
    fprintf(fid, '\nAll bones - %d/%d (%8.6g%%) non-HSCs in significant clones belong to HSC-derived clones (p=%3.2e)\n', ...
        non_HSC_in_clone, non_HSC, non_HSC_in_clone/non_HSC*100, EB.p_derived_in_hsc_derived);
         
    fprintf(fid, '\nNon-uniformity in HSC-rooted clone size: p=%3.2e\n', EB.csd.pval);
    
    bone_membership_of_sig_clones = sum(logical(combined.allele_breakdown_by_sample(combined.sig_alleles.all,:)),2);
    sig_clone_in_all_bones = combined.sig_alleles.all(bone_membership_of_sig_clones==4);
    fprintf(fid, '\nNumber of significant CARLIN alleles that show up in all bones: %d\n', length(sig_clone_in_all_bones));
    fprintf(fid, '\n%% edited cells in significant clones found in all bones: %6.5g\n', ...
        sum(sum(combined.allele_breakdown_by_sample(sig_clone_in_all_bones,:),2))/sum(sum(combined.allele_breakdown_by_sample(2:end,:),2))*100);
    
    highlight_clone = combined.sig_alleles.all(find(bone_membership_of_sig_clones==2,1,'first'));
    [~, highlight_clone_bone_name, highlight_clone_bone_vals] = find(combined.allele_breakdown_by_sample(highlight_clone,:));
    fprintf(fid, '\nClone %d appears in %d cells in %s and %d cells in %s (p=%3.2e)\n', highlight_clone, ...
        highlight_clone_bone_vals(1), sample_names{highlight_clone_bone_name(1)}, ...
        highlight_clone_bone_vals(2), sample_names{highlight_clone_bone_name(2)}, ...
        bb_pval(highlight_clone)/bonferroni_correction_factor);
    
    highlight_clone = combined.sig_alleles.all(find(bone_membership_of_sig_clones==4,1,'first'));
    [~, highlight_clone_bone_name, highlight_clone_bone_vals] = find(combined.allele_breakdown_by_sample(highlight_clone,:));
    fprintf(fid, '\nClone %d appears in [%d,%d,%d,%d] cells in [%s,%s,%s,%s] (p=%3.2e)\n', highlight_clone, ...
        highlight_clone_bone_vals(1), highlight_clone_bone_vals(2), highlight_clone_bone_vals(3), highlight_clone_bone_vals(4), ...
         sample_names{highlight_clone_bone_name(1)}, sample_names{highlight_clone_bone_name(2)}, ...
         sample_names{highlight_clone_bone_name(3)}, sample_names{highlight_clone_bone_name(4)}, ...
         bb_pval(highlight_clone)/bonferroni_correction_factor);

     
    %% 5FU
    
    fprintf(fid, '\n5FU\n');

    load(sprintf('%s/5FU/Analysis.mat', results_path));
    N_cells = cellfun(@(x) x.summary.N.reference_tags, samples);    
    fprintf(fid, '\nRange cells for 5FU: [%d, %d]\n', min(N_cells), max(N_cells));
    
    fprintf(fid, '\nTotal number of Louvain clusters: %d\n', max(pooled.louvain)+1);
    
    summaries = cellfun(@(x) x.summary, samples, 'un', false);
    [~, sample_map, ~] = ExperimentSummary.FromMerge(vertcat(summaries{:}));
    rare_alleles = cellfun(@(x,y) x(y.sig_alleles.all), sample_map, samples, 'un', false);
    rare_alleles = vertcat(rare_alleles{:});
    
    fprintf(fid, '\nAcross 5 mice: %d/%d (%8.6g%%) rare alleles found in only one mouse\n', ...
        length(unique(rare_alleles)), length(rare_alleles), length(unique(rare_alleles))/length(rare_alleles)*100);

    load(sprintf('%s/SC/Stats.mat', results_path), 'FU5', 'ctrl');
    
    fprintf(fid, '\nReduction in number of clones: p=%3.2e\n', FU5.nc{1,1}.pval);
    
    fprintf(fid, '\nControl: %d/%d (%8.6g%%) significant alleles are parent clones representing %d/%d (%8.6g%%) cells marked with significant alleles\n', ...
        ctrl.clones.hsc_derived, ctrl.clones.sig, ctrl.clones.hsc_derived/ctrl.clones.sig*100, ...
        ctrl.cells.in_hsc_derived, ctrl.cells.in_sig, ctrl.cells.in_hsc_derived/ctrl.cells.in_sig*100);
        
    fprintf(fid, '\n5FU: %d/%d (%8.6g%%) significant alleles are parent clones representing %d/%d (%8.6g%%) cells marked with significant alleles\n', ...
        FU5.clones.hsc_derived, FU5.clones.sig, FU5.clones.hsc_derived/FU5.clones.sig*100, ...
        FU5.cells.in_hsc_derived, FU5.cells.in_sig, FU5.cells.in_hsc_derived/FU5.cells.in_sig*100);
    
    fprintf(fid, '\np-Value increased fraction of parent clones: %3.2e\n', FU5.p_clone_is_hsc_derived);
    
    fprintf(fid, '\np-Value increased fraction of cells in parent clones: %3.2e\n', FU5.p_cell_in_hsc_derived);
        
    fprintf(fid, '\nIncreased HSC-rooted clone size in 5FU: p=%3.2e\n', FU5.acs{1,1}.pval);
    
    hsc_rooted_cs = cellfun(@(x) x.summary.allele_freqs(union(x.sig_alleles.hsc_derived_allele, x.sig_alleles.hsc_only_allele)), samples(1:2), 'un', false);
        
    top_cells = sum(hsc_rooted_cs{pos_samples(1)}(1:12));
    all_cells = sum(hsc_rooted_cs{pos_samples(1)});
    fprintf(fid, '\n5FU Mouse 1: Top 12 of %d HSC-rooted clones account for %d/%d (%8.6g%%) cells in HSC-rooted clones (p=%3.2e)\n', ...
        length(hsc_rooted_cs{pos_samples(1)}), top_cells, all_cells, top_cells/all_cells*100, FU5.csd{1}.pval);
    
    top_cells = sum(hsc_rooted_cs{pos_samples(2)}(1:4));
    all_cells = sum(hsc_rooted_cs{pos_samples(2)});
    fprintf(fid, '\n5FU Mouse 2: Top 4 of %d HSC-rooted clones account for %d/%d (%8.6g%%) cells in HSC-rooted clones (p=%3.2e)\n', ...
        length(hsc_rooted_cs{pos_samples(2)}), top_cells, all_cells, top_cells/all_cells*100, FU5.csd{2}.pval);

    load(sprintf('%s/SC/Stats.mat', results_path), 'FU5_ctrl');
    
    fprintf(fid, '\nNumber of [parent/childless] HSCs cells: %d/%d\n', FU5_ctrl.cells.hsc_parent, FU5_ctrl.cells.hsc_childless);
    
    fprintf(fid, '\nParent HSCs over-represented in parent cluster: [z=%.2f, p=%3.2e]\n', ...
        FU5_ctrl.z_parent_HSC_in_parent_cluster, FU5_ctrl.p_parent_HSC_in_parent_cluster);
    
    fprintf(fid, '\nProliferating HSCs over-represented in parent cluster: [z=%.2f, p=%3.2e]\n', ...
        FU5_ctrl.z_high_prolif_in_parent_cluster, FU5_ctrl.p_high_prolif_in_parent_cluster);
    
    fprintf(fid, '\nProliferating HSCs over-represented in parent HSCs: [z=%.2f, p=%3.2e]\n', ...
        FU5_ctrl.z_high_prolif_in_parent_HSC, FU5_ctrl.p_high_prolif_in_parent_HSC);
   
    
    %% Supplemental
    
    fprintf(fid, '\nSUPPLEMENTAL\n');
    bank = pos_bank.bank;
    
    p = [bank.model.extensive.rates(2:end); ones(bank.model.alleles.unobserved,1)*bank.model.extensive.mean_unobs_rate];
    p = p/sum(p);
    H = -sum(p.*log2(p));
    fprintf(fid, '\nEffective alleles in bank: %d\n', 2^H);
    
    max_rate = Bank.max_rate_for_sig_level(0.05, 5000);
    max_rate = max_rate{1};
    fprintf(fid, '\nFor 5000 induced cells, %3.2g%% of alleles are 0.05 significant representing %3.2g%% of cells\n', ...
        (sum(bank.model.intensive.rates < max_rate)+bank.model.alleles.unobserved)/bank.model.alleles.estimated*100, ...
        (sum(bank.model.intensive.rates(bank.model.intensive.rates < max_rate))+bank.model.intensive.mean_unobs_rate*bank.model.alleles.unobserved)/...
        (sum(bank.model.intensive.rates(2:end))+bank.model.intensive.mean_unobs_rate*bank.model.alleles.unobserved)*100);
    
    fclose(fid);
    
end