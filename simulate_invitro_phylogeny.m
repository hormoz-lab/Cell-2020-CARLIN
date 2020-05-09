function simulate_invitro_phylogeny(analysis_dir, results_dir)

    outdir = sprintf('%s/Trees/InVitroPhylogeny', results_dir);

    for i = 1:8
        samples{i} = load(sprintf('%s/InVitroPhylogeny/Gen2/%d/Summary.mat', analysis_dir, i), 'summary');
        samples{i} = samples{i}.summary;
    end
    samples = samples';
    filter_params = struct('epval', 0, 'cpval', 0, 'pos', false, 'neg', false, 'ambiguous', false, 'potential', 0);
    rng(1973);
    
    tissue_reconstruction(samples, outdir, filter_params);

    load(sprintf('%s/Lineage.mat', outdir));

    % 1. Load Sanger alleles

    sanger = cell(N_samples,1);
    for i = 1:N_samples
        sanger{i} = load(sprintf('Analyzed/InVitroPhylogeny/Gen1/%d/Summary.mat', i), 'summary');
        sanger{i} = sanger{i}.summary.alleles{1};
    end

    % 2. Generate table of mutations found in Sanger colonies

    sanger_events = cellfun(@Mutation.identify_Cas9_events, sanger, 'un', false);
    N_sanger_muts = cellfun(@length, sanger_events);
    sanger_events = vertcat(sanger_events{:});
    sanger_events = arrayfun(@(i) sanger_events(i).annotate(false), [1:length(sanger_events)]', 'un', false);
    [master_sanger_events, ~, which_event] = unique(sanger_events, 'stable');
    master_event_table = false(length(master_sanger_events),N_samples);
    master_event_table(sub2ind(size(master_event_table), which_event, repelem([1:N_samples]', N_sanger_muts))) = true;
    master_event_table = master_event_table.*[1:N_samples];

    %3. Locate the Gen1 Sanger sequences in the Gen2 bulk sequences.

    sanger_allele = NaN(N_samples,1);

    [is, where] = ismember(cellfun(@(x) degap(x.get_seq), sanger, 'un', false), ...
                           cellfun(@(x) degap(x.get_seq), aggregated_alleles, 'un', false));

    sanger_allele(is) = where(is);
    sanger_allele = [1; sanger_allele];

    % 4. Find Sanger progenitor distribution for each allele across tree ensemble.
    %    We define the Sanger progenitor of an allele to be the first internal
    %    node on the path from leaf back to the root, which corresponds to a
    %    Sanger allele.

    % Stores ancestors from leaf to root
    retrace_progenitor = cellfun(@(x,y) fliplr(x(y)), allele_at_node, nodes_on_path_to_leaf, 'un', false);
    
    [is_sanger_progenitor, which_sanger_progenitor] = cellfun(@(x) arrayfun(@(i) ismember(x(i,:)', sanger_allele), ...
                                                                   [1:N_alleles], 'un', false), retrace_progenitor, 'un', false);
                                                               
    is_sanger_progenitor = vertcat(is_sanger_progenitor{:});
    which_sanger_progenitor = vertcat(which_sanger_progenitor{:});
    sanger_progenitor = cellfun(@(x,y) y(find(x, 1, 'first')), is_sanger_progenitor, which_sanger_progenitor)';    
    sanger_progenitor_dist = accumarray([repmat([1:N_alleles]', [N_sim,1]), sanger_progenitor(:)], 1, [N_alleles, N_samples+1]);
        
    tn = arrayfun(@(i) (sanger_progenitor(1,i)==1)*allele_freq(1)/sum(allele_freq), [1:N_sim]');
    fn = arrayfun(@(i) sum((sanger_progenitor(2:end,i)==1).*allele_freq(2:end))/sum(allele_freq), [1:N_sim]');
    tp = arrayfun(@(i) sum(sum(full(sparse([1:N_alleles-1]', sanger_progenitor(2:end,i), sanger_progenitor(2:end,i)>1, N_alleles-1, N_samples+1)...
                            .*[zeros(N_alleles-1,1) allele_breakdown_by_sample(2:end,:)])))/sum(allele_freq), [1:N_sim]');
    fp = 1-tn-fn-tp;
    
    save(sprintf('%s/Analysis.mat', outdir), ...
        'sanger', 'sanger_events', 'N_sanger_muts', 'master_sanger_events', 'which_event', 'master_event_table', ...
        'sanger_allele', 'sanger_progenitor', 'sanger_progenitor_dist', 'fp', 'fn', 'tp', 'tn');
end
