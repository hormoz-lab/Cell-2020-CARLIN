function [x, y, desired_order] = lineage_tree_layout(adj_list, allele_at_node, allele_pop, reorder_hook, hook_priority)

    N_alleles = find(adj_list==0)-1;

    assert(length(allele_pop) == N_alleles && length(reorder_hook) == N_alleles);
    assert(all(allele_pop>0));
    assert(adj_list(N_alleles+1)==0 && allele_at_node(N_alleles+1)==1);    
    assert(all(ismember(reorder_hook, [1:length(hook_priority)])));
    assert(length(nonzeros(hook_priority))==length(unique(nonzeros(hook_priority))));
        
    [node_depth, node_height, chain, clade_pop] = compute_node_properties(adj_list, allele_at_node, allele_pop);

    max_depth = max(node_depth)+1;
    
    % Ordering rules
    %   a. First according to reorder hook priority applied to each leaf's
    %      progenitor according to depth
    %   b. Second according to leaf depth
    %   c. Third according to common parentage
    %   d. Last according to population
    
    % allele_at_node may not necessarily be contiguously numbered if the
    % adj list being passed is from a pruned tree    
    [~, allele_at_node] = ismember(allele_at_node, allele_at_node(1:N_alleles));
    
    priority_vector = zeros(N_alleles, 4*max_depth);
    priority_vector(:,1:4:end) = hook_priority(reorder_hook(allele_at_node(chain)));
    priority_vector(:,2:4:end)   = node_depth(chain);
    priority_vector(:,3:4:end)   = chain;
    priority_vector(:,4:4:end)   = clade_pop(chain);
    
    [~, desired_order] = sortrows(priority_vector, 'descend');

    x = node_height;
    y = NaN(size(x));
    y(desired_order) = [1:N_alleles]';
    
    for h = 1:max(node_height)
        for n = find(node_height==h)'        
            y(n) = mean(y(chain(:,node_depth(n)+1)==n));
        end
    end
        
    assert(~any(isnan(y)));
    
    desired_order = [desired_order; [N_alleles+1:length(adj_list)]'];
    
end