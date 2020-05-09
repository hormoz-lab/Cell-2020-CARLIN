function make_algorithm_subplots(results_dir)

    outdir = sprintf('%s/Figures/Algorithm', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    load(sprintf('%s/Banks/Protocol2/RNA/PosDox/Bank.mat', results_dir));
    
    ref = CARLIN_def.getInstance;
    
    Primer5 = ref.seq.Primer5;
    Primer3 = ref.seq.Primer3;
    
    % Flank with primers. Other papers do this to avoid ragged ends,
    % and we want to make a fair comparison.
    [~, nw_aligned] = nwalign(cellfun(@(x) [Primer5 degap(x.get_seq) Primer3], bank.summary.alleles, 'un', false), ...
                              [Primer5 ref.seq.CARLIN Primer3], 'Alphabet', 'NT', 'GapOpen', 10, 'ExtendGap', 0.5);
    [ref_inds, ref_inde] = cellfun(@(x) deal(find(x(3,:)~='-', ref.width.Primer5+1, 'first'), ...
                                             find(x(3,:)~='-', ref.width.Primer3+1, 'last')), nw_aligned, 'un', false);
    ref_inds = cellfun(@(x) x(end), ref_inds);
    ref_inde = cellfun(@(x) x(1), ref_inde);    
    nw_aligned = arrayfun(@(i) nw_aligned{i}(:,ref_inds(i):ref_inde(i)), [1:size(nw_aligned,1)]', 'un', false);
    nw_aligned = cellfun(@(x) CARLIN_def.desemble_sequence(x(1,:), x(3,:)), nw_aligned, 'un', false);    
    
    cas9_aligned = bank.summary.alleles;
    N_trunc = 10;
    
	nw_split   = show_alignment_frame(cellfun(@Mutation.identify_sequence_events, nw_aligned, 'un', false), ref, N_trunc);
    cas9_split   = show_alignment_frame(cellfun(@Mutation.identify_sequence_events, cas9_aligned, 'un', false), ref, N_trunc);
    
	nw_grouped = show_alignment_frame(cellfun(@Mutation.identify_Cas9_events, nw_aligned, 'un', false), ref, N_trunc);
    cas9_grouped = show_alignment_frame(cellfun(@Mutation.identify_Cas9_events, cas9_aligned, 'un', false), ref, N_trunc);
    
    nw_sc_cas9_aligned = cellfun(@(x) CARLIN_def.nwalign_score(x), cas9_aligned);
    nw_sc_nw_aligned   = cellfun(@(x) CARLIN_def.nwalign_score(x), nw_aligned);
    nw_sc_diff_aligned = nw_sc_nw_aligned-nw_sc_cas9_aligned;
    assert(all(nw_sc_diff_aligned>=0));
    [nw_sc_diff, nw_sc_diff_CDF] = unique(sort(nw_sc_diff_aligned), 'last');
    plot_alignment_score_difference(nw_sc_diff, nw_sc_diff_CDF/length(nw_sc_diff_aligned));
    paper_print(sprintf('%s/%s', outdir, 'ScoreComparison'));
    
    series_names = {'NW Unmerged'; sprintf('NW Merged'); 'CARLIN Unmerged'; sprintf('CARLIN Merged')};

    plot_mutation_alignment_count([nw_split.N.indel nw_grouped.N.indel cas9_split.N.indel cas9_grouped.N.indel], ...
                                  'Indel Frequency', series_names, true);
    paper_print(sprintf('%s/%s', outdir, 'IndelFrequency'));

    plot_mutation_alignment_count([nw_split.N.del   nw_grouped.N.del   cas9_split.N.del   cas9_grouped.N.del], ...
                                  'Deletion Frequency', series_names);
    paper_print(sprintf('%s/%s', outdir, 'DeletionFrequency'));    
    
    plot_mutation_alignment_count([nw_split.N.ins   nw_grouped.N.ins   cas9_split.N.ins   cas9_grouped.N.ins], ...
                                  'Insertion Frequency', series_names);
    paper_print(sprintf('%s/%s', outdir, 'InsertionFrequency'));

    plot_mutation_alignment_frame([nw_split.start_bp.indel nw_grouped.start_bp.indel cas9_split.start_bp.indel cas9_grouped.start_bp.indel], ...
                                  'Indel Starting BP', series_names);
    paper_print(sprintf('%s/%s', outdir, 'IndelStart'));
                              
    plot_mutation_alignment_frame([nw_split.end_bp.indel nw_grouped.end_bp.indel cas9_split.end_bp.indel cas9_grouped.end_bp.indel], ...
                                  'Indel Ending BP', series_names);
    paper_print(sprintf('%s/%s', outdir, 'IndelEnd'));
    
    plot_mutation_alignment_frame([nw_split.start_bp.del nw_grouped.start_bp.del cas9_split.start_bp.del cas9_grouped.start_bp.del], ...
                                  'Deletion Starting BP', series_names);
    paper_print(sprintf('%s/%s', outdir, 'DeletionStart'));
                              
    plot_mutation_alignment_frame([nw_split.end_bp.del nw_grouped.end_bp.del cas9_split.end_bp.del cas9_grouped.end_bp.del], ...
                                  'Deletion Ending BP', series_names);
    paper_print(sprintf('%s/%s', outdir, 'DeletionEnd'));
    
    plot_mutation_alignment_frame([nw_split.start_bp.ins nw_grouped.start_bp.ins cas9_split.start_bp.ins cas9_grouped.start_bp.ins], ...
                                  'Insertion Site BP', series_names);
    paper_print(sprintf('%s/%s', outdir, 'InsertionSite'));
    
    close all;
    
end

function out = show_alignment_frame(mut_list, ref, N_trunc)
    
    L_mut = cellfun(@length, mut_list);
    
    mut_list = vertcat(mut_list{:});
    
    start_loc = arrayfun(@(i) CARLIN_def.locate(ref, mut_list(i).loc_start, ref.bounds.ordered), [1:length(mut_list)]');
    end_loc   = arrayfun(@(i) CARLIN_def.locate(ref, mut_list(i).loc_end,   ref.bounds.ordered), [1:length(mut_list)]');    
    
    % Restrict analysis to mutations that don't get into prefix and
    % postfix, which fall outside alignment frame.
        
    mut_type  = vertcat(mut_list.type);
    
    out.N.del   = cellfun(@sum, mat2cell(mut_type == 'D', L_mut));
    out.N.ins   = cellfun(@sum, mat2cell(mut_type == 'I', L_mut));
    out.N.indel = cellfun(@sum, mat2cell(mut_type == 'C' | mut_type == 'M', L_mut));
    
    out.N.del = accumarray(out.N.del+1, 1);
    out.N.del(N_trunc+1) = sum(out.N.del(N_trunc+1:end));
    out.N.del = out.N.del(1:N_trunc+1);
    
    out.N.ins = accumarray(out.N.ins+1, 1);
    out.N.ins(N_trunc+1) = sum(out.N.ins(N_trunc+1:end));
    out.N.ins = out.N.ins(1:N_trunc+1);
    
    out.N.indel = accumarray(out.N.indel+1, 1);
    out.N.indel(N_trunc+1) = sum(out.N.indel(N_trunc+1:end));
    out.N.indel = out.N.indel(1:N_trunc+1);
    
    valid_sites = ~ (ismember({start_loc.type}', {'prefix'; 'postfix'}) | ismember({end_loc.type}', {'prefix'; 'postfix'}));
    start_loc = start_loc(valid_sites);
    end_loc   = end_loc(valid_sites);   
    mut_type  = mut_type(valid_sites);
     
    start_bp = vertcat(start_loc.pos);
    [is, where] = ismember({start_loc.type}', {'consites'; 'cutsites'; 'pams'});
    assert(all(is));
    start_bp(where==2) = vertcat(start_loc(where==2).pos) + ref.width.consite;
    start_bp(where==3) = vertcat(start_loc(where==3).pos) + ref.width.segment;
    
    end_bp = vertcat(end_loc.pos);
    [is, where] = ismember({end_loc.type}, {'consites'; 'cutsites'; 'pams'});
    assert(all(is));
    end_bp(where==2) = vertcat(end_loc(where==2).pos) + ref.width.consite;
    end_bp(where==3) = vertcat(end_loc(where==3).pos) + ref.width.segment;
    
    W = ref.width.segment + ref.width.pam;
    
    out.start_bp.del   = accumarray(start_bp(mut_type == 'D'), 1, [W, 1]);
    out.start_bp.ins   = accumarray(start_bp(mut_type == 'I'), 1, [W, 1]);
    out.start_bp.indel = accumarray(start_bp(mut_type == 'C' | mut_type == 'M'), 1, [W, 1]);
    
    out.end_bp.del     = accumarray(end_bp(mut_type == 'D'), 1, [W, 1]);
    out.end_bp.ins     = accumarray(end_bp(mut_type == 'I'), 1, [W, 1]);
    out.end_bp.indel   = accumarray(end_bp(mut_type == 'C' | mut_type == 'M'), 1, [W, 1]);
    
    assert(isequal(out.start_bp.ins, out.end_bp.ins));
    out.end_bp.ins = zeros(size(out.end_bp.ins));
    
end

 
