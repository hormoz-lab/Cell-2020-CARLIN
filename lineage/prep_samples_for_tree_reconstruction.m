function [raw_combined, allele_breakdown_by_sample, aggregated_alleles] = prep_samples_for_tree_reconstruction(samples, filter_params)

    % Prune alleles according to the parameters outlined in filter_params

    N_samples = length(samples);
    [raw_combined.summary, raw_combined.sample_map, raw_combined.allele_breakdown_by_sample] = ExperimentSummary.FromMerge(vertcat(samples{:}));    
        
    % Keep the template sequence first. Subsequent code relies on it being
    % the root of the tree.
    
    [allele_breakdown_by_sample, reorder] = sortrows(raw_combined.allele_breakdown_by_sample, [1:N_samples]', 'descend');    
    aggregated_alleles = raw_combined.summary.alleles(reorder);    
    template_ind = find(strcmp(cellfun(@(x) degap(x.get_seq), aggregated_alleles, 'un', false), CARLIN_def.getInstance.seq.CARLIN));
    assert(~isempty(template_ind), 'Tissue reconstruction needs template in allele pool');
    
    allele_breakdown_by_sample = [allele_breakdown_by_sample(template_ind,:); ...
                                  allele_breakdown_by_sample(1:template_ind-1,:); ...
                                  allele_breakdown_by_sample(template_ind+1:end,:)];
                              
    aggregated_alleles = [aggregated_alleles(template_ind); aggregated_alleles(1:template_ind-1); aggregated_alleles(template_ind+1:end)];
    
    if (filter_params.pos)
        bank1 = load('Results/Banks/Protocol1/RNA/PosDox/Bank.mat');
        bank1 = bank1.bank;
        is1 = ismember(cellfun(@(x) degap(x.get_seq), aggregated_alleles, 'un', false), ...
                          cellfun(@(x) degap(x.get_seq), bank1.summary.alleles, 'un', false));                      
        bank2 = load('Results/Banks/Protocol2/RNA/PosDox/Bank.mat');
        bank2 = bank2.bank;        
        is2 = ismember(cellfun(@(x) degap(x.get_seq), aggregated_alleles, 'un', false), ...
                       cellfun(@(x) degap(x.get_seq), bank2.summary.alleles, 'un', false));        
        rare = [1; find(~is1 & ~is2)];
        allele_breakdown_by_sample = allele_breakdown_by_sample(rare,:);
        aggregated_alleles = aggregated_alleles(rare,:);    
    end
    
    if (filter_params.cpval > 0)
        load('Results/Banks/Protocol2/RNA/PosDox/Bank.mat');
        p = bank.compute_clonal_pvalue(aggregated_alleles(2:end), sum(sum(allele_breakdown_by_sample(2:end,:),2)));
        ind = [1; find(p <= filter_params.cpval)+1];
        allele_breakdown_by_sample = allele_breakdown_by_sample(ind,:);
        aggregated_alleles = aggregated_alleles(ind);
    end
    
    if (filter_params.epval > 0)
        bank1 = load('Results/Banks/Protocol1/RNA/PosDox/Bank.mat');
        bank1 = bank1.bank;
        p1 = bank1.compute_frequency_pvalue(aggregated_alleles(2:end), sum(allele_breakdown_by_sample(2:end,:),2));
        bank2 = load('Results/Banks/Protocol2/RNA/PosDox/Bank.mat');
        bank2 = bank2.bank;
        p2 = bank2.compute_frequency_pvalue(aggregated_alleles(2:end), sum(allele_breakdown_by_sample(2:end,:),2));
        p = max(p1, p2);
        pcrit = benjamini_hochberg(p, filter_params.epval);
        ind = [1; find(p <= pcrit)+1];
        allele_breakdown_by_sample = allele_breakdown_by_sample(ind,:);
        aggregated_alleles = aggregated_alleles(ind);
    end
    
    if (filter_params.neg)
        bank1 = load('Results/Banks/Protocol1/RNA/NegDox/Bank.mat');
        bank1 = bank1.bank;
        is1 = ismember(cellfun(@(x) degap(x.get_seq), aggregated_alleles, 'un', false), ...
                          cellfun(@(x) degap(x.get_seq), bank1.summary.alleles, 'un', false));                      
        bank2 = load('Results/Banks/Protocol2/RNA/NegDox/Bank.mat');
        bank2 = bank2.bank;        
        is2 = ismember(cellfun(@(x) degap(x.get_seq), aggregated_alleles, 'un', false), ...
                       cellfun(@(x) degap(x.get_seq), bank2.summary.alleles, 'un', false));        
        rare = [1; find(~is1 & ~is2)];
        allele_breakdown_by_sample = allele_breakdown_by_sample(rare,:);
        aggregated_alleles = aggregated_alleles(rare,:);    
    end
    
    if (filter_params.ambiguous)
        muts = cellfun(@Mutation.identify_Cas9_events, aggregated_alleles, 'un', false);
        ambiguous = cellfun(@allele_has_ambiguous_deletions, muts);        
        allele_breakdown_by_sample = allele_breakdown_by_sample(~ambiguous,:);
        aggregated_alleles = aggregated_alleles(~ambiguous,:);        
    end
    
    if (filter_params.potential > 0)
        CP = cellfun(@(x) CARLIN_def.getInstance.N.segments-length(Mutation.find_modified_sites(x)), aggregated_alleles);
        aggregated_alleles = aggregated_alleles(CP>=cull_potential);
        allele_breakdown_by_sample = allele_breakdown_by_sample(CP>=cull_potential,:);
    end
    
    assert(isequal(aggregated_alleles{1}.get_seq, CARLIN_def.getInstance.seq.CARLIN));
    
end