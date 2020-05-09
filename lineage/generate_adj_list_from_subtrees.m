function [adj_list, allele_at_node, N_pruned_alleles] = ...
    generate_adj_list_from_subtrees(stable_depthX_subtrees, depthX_subtree_clade_pop, which_stat, threshold)

    % 1. Gen1 alleles need to be stable and exceed a stat threshold, or be
    %    the parent of stable Gen2 subtree.

    gen1_al = unique([stable_depthX_subtrees{1}(depthX_subtree_clade_pop{1}.(which_stat)>=threshold); stable_depthX_subtrees{2}(:,1)]);
    
    % 2. Initialize child leaves and set initial parent to be template
    %    root located at (N_stable_alleles+1).
    
    allele_at_node = [gen1_al; stable_depthX_subtrees{2}(:,2); stable_depthX_subtrees{3}(:,3); 1];
    N_pruned_alleles = length(allele_at_node)-1;
    adj_list = [(N_pruned_alleles+1)*ones(N_pruned_alleles,1); 0];
    
    % 3. Create internal nodes for leaf alleles that are also Gen2
    % progenitors
    
    [~, where_gen2_in_gen1] = ismember(stable_depthX_subtrees{2}(:,1), allele_at_node);    
    [where_gen2_in_gen1, ~, which_parent] = unique(where_gen2_in_gen1);
    gen1_new_nodes = length(adj_list)+[1:length(where_gen2_in_gen1)]';
    gen1_old_head = adj_list(where_gen2_in_gen1);
    adj_list(where_gen2_in_gen1) = gen1_new_nodes;
    adj_list(length(gen1_al)+[1:length(stable_depthX_subtrees{2}(:,2))]') = gen1_new_nodes(which_parent);
    adj_list = [adj_list; gen1_old_head];
    allele_at_node = [allele_at_node; allele_at_node(where_gen2_in_gen1)];
    
    % 4. Create internal nodes for Gen2 alleles that are also Gen3
    % progenitors

    [~, where_gen3_in_gen2] = ismember(stable_depthX_subtrees{3}(:,2), allele_at_node);
    [where_gen3_in_gen2, ~, which_parent] = unique(where_gen3_in_gen2);
    gen2_new_nodes = length(adj_list)+[1:length(where_gen3_in_gen2)]';
    gen2_old_head = adj_list(where_gen3_in_gen2);
    adj_list(where_gen3_in_gen2) = gen2_new_nodes;
    adj_list(length(gen1_al)+length(stable_depthX_subtrees{2}(:,2))+[1:length(stable_depthX_subtrees{3}(:,3))]') = gen2_new_nodes(which_parent);
    adj_list = [adj_list; gen2_old_head];
    allele_at_node = [allele_at_node; allele_at_node(where_gen3_in_gen2)];

end