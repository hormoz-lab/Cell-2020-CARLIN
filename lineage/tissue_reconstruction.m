function tissue_reconstruction(samples, outpath, results_path, filter_params, N_sim)

    if (nargin < 5)
        N_sim = 10000;
    end
  
    % 1. Filter alleles according to filter_params and merge them across samples
    
    [raw_combined, allele_breakdown_by_sample, aggregated_alleles] = prep_samples_for_tree_reconstruction(samples, filter_params, results_path);
    [N_alleles, N_samples] = size(allele_breakdown_by_sample);
    
    % 2. Compute pairwise lineage relationships between all alleles
    
    [mut_matrix, relationship, mutations, type, preserved, mut_events] = compute_allele_hierarchy(aggregated_alleles);
    
    % 3. Create a binary adjacency matrix based on pairwise information
    
    allele_freq = sum(allele_breakdown_by_sample,2);
    adj = set_adjacency_matrix(mut_matrix, relationship, mutations, type, preserved, allele_freq);
    
    % 4. Do probabilistic simulations
    
    adj_list = cell(N_sim,1);
    allele_at_node = cell(N_sim,1);
    parfor i = 1:N_sim
        [i N_sim]
        [adj_list{i}, allele_at_node{i}] = sim_tree([N_alleles+1; NaN(N_alleles-1,1); 0], [1; NaN(N_alleles-1,1); 1], adj, 1);
    end
    
    % 5. Compute the height, depth, path and clade population of each node
    
    fprintf('Computing node properties\n');
    [node_depth, node_height, nodes_on_path_to_leaf, node_subclade_pop] = ...
        cellfun(@(x,y) compute_node_properties(x, y, allele_freq), adj_list, allele_at_node, 'un', false);
   
    % 6. Compute tree isomorphisms and redundancies
    
    isomorphic = check_tree_isomorphisms_allele(adj, adj_list, allele_at_node);
    [redundancy, unique_trees] = sort(accumarray(isomorphic, 1), 'descend');
   
    % 7. Compute stable tree
    
    [stable_adj_list, stable_allele_at_node, N_stable_alleles, stable_depthX_subtrees, N_omnipresent_depthX_subtrees, ...
     depthX_subtree, depthX_subtree_instances, depthX_subtree_pop, depthX_subtree_clade_pop, is_allele_at_depthX, depthX_alleles] = ...          
          compute_stable_tree(N_sim, N_alleles, N_samples, nodes_on_path_to_leaf, allele_at_node, allele_breakdown_by_sample);
  
    if (~exist(outpath, 'dir'))
        mkdir(outpath);
    end
    
    save(sprintf('%s/Lineage.mat', outpath));
    
end