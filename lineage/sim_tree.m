function [adj_list, node_allele] = sim_tree(adj_list, node_allele, adj_matrix, current_node)

    N_alleles = size(adj_matrix,1);
    assert(~any(adj_list(1:N_alleles)==current_node));
    
    unplaced_alleles = isnan(adj_list(1:N_alleles));
    assert(~unplaced_alleles(current_node));    
    
    alleles_to_consider = find(unplaced_alleles & adj_matrix(current_node,:)');
    
    while (~isempty(alleles_to_consider))
        num_children = sum(adj_matrix(alleles_to_consider,:) & unplaced_alleles',2);
        assert(all(num_children));
        allele_freq = cumsum(num_children)/sum(num_children);
        p = rand;
        allele_chosen = find(allele_freq >= p, 1, 'first');
        if (num_children(allele_chosen)==1)
            adj_list(alleles_to_consider(allele_chosen)) = adj_list(current_node);
            node_allele(alleles_to_consider(allele_chosen)) = alleles_to_consider(allele_chosen);
        else
            adj_list = [adj_list; adj_list(current_node)];
            node_allele = [node_allele; alleles_to_consider(allele_chosen)];
            adj_list(alleles_to_consider(allele_chosen)) = size(adj_list,1);
            node_allele(alleles_to_consider(allele_chosen)) = alleles_to_consider(allele_chosen);
            [adj_list, node_allele] = sim_tree(adj_list, node_allele, adj_matrix, alleles_to_consider(allele_chosen));            
        end
        unplaced_alleles = isnan(adj_list(1:N_alleles));
        alleles_to_consider(ismembc(alleles_to_consider, find(~unplaced_alleles))) = [];
    end
    return;
end
