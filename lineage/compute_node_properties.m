function [node_depth, node_height, nodes_on_path_to_leaf, node_subclade_pop] = compute_node_properties(adj_list, allele_at_node, allele_pop)

    N_alleles = find(adj_list==0)-1;

    assert(size(allele_pop,1) == N_alleles);
    N_samples = size(allele_pop,2);
    assert(all(sum(allele_pop,2)>0));
    assert(adj_list(N_alleles+1)==0 && allele_at_node(N_alleles+1)==1);    
    
    parent_list = 0;
    left_to_classify = 1:length(adj_list);    
    node_depth = NaN(size(adj_list));
    
    depthX = 0;
    while (~isempty(left_to_classify))        
        is_depthX = left_to_classify(ismember(adj_list(left_to_classify), parent_list));
        node_depth(is_depthX) = depthX;
        left_to_classify = setdiff(left_to_classify, is_depthX);
        parent_list = is_depthX;
        depthX = depthX + 1;
    end
    
    assert(~any(isnan(node_depth)));
    
    is_heightX = 1:N_alleles;
    node_height = NaN(size(adj_list));
    
    heightX = 0;
    while (~isempty(is_heightX))
        node_height(is_heightX) = heightX;
        is_heightX = setdiff(unique(adj_list(is_heightX)), 0);
        heightX = heightX + 1;
    end
        
    assert(~any(isnan(node_height)));    
    
    max_depth = max(max(node_depth),3);
    
    nodes_on_path_to_leaf = NaN(N_alleles,max_depth+1);
    node_subclade_pop = zeros(size(allele_at_node,1), N_samples);
    
    for k = 1:N_alleles            
        cur = k;
        ind = 1;
        nodes_on_path_to_leaf(k,ind) = cur;
        node_subclade_pop(cur,:) = allele_pop(k,:);
        while (adj_list(cur) ~= 0)
            cur = adj_list(cur);
            ind = ind+1;
            nodes_on_path_to_leaf(k,ind) = cur;            
            node_subclade_pop(cur,:) = node_subclade_pop(cur,:)+allele_pop(k,:);
        end
        nodes_on_path_to_leaf(k,1:ind) = fliplr(nodes_on_path_to_leaf(k,1:ind));
        nodes_on_path_to_leaf(k,ind+1:end) = nodes_on_path_to_leaf(k,ind);
    end
    
    assert(~any(isnan(nodes_on_path_to_leaf(:))));        
    
    assert(isequal(node_subclade_pop(N_alleles+1,:), sum(allele_pop,1)));
    
end