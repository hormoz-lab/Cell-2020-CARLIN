function isomorphic = check_tree_isomorphisms_allele(adj_true, adj_tree, allele_at_node)

    fprintf('Checking tree isomorphisms\n');

    N_alleles = size(adj_true,1);
    N_sim = size(adj_tree,1);
    
    % One node is marked as root by having its parent be 0. This is always
    % the template, so readjust this so we can non-zero index.    
    adj_tree = cellfun(@(x) max(x, 1), adj_tree, 'un', false);
    
    % Make adjacency matrix from adjacency list
    adj_tree = arrayfun(@(i) sparse(logical(accumarray([allele_at_node{i}(adj_tree{i}), allele_at_node{i}], true, ...
                                    [N_alleles, N_alleles]))), [1:N_sim]', 'un', false);
                                
    % Check if adjacency matrices are unique. Reshape to [N_sim,
    % N_allele*N_allele] matrix so we can use MATLAB's unique(...'rows')
    % function
    
    adj_tree = horzcat(adj_tree{:});        % [Adj_1, Adj_2, ..., Adj_Nsim] 
    
    % [N_sim, N_allele x N_allele]
    adj_tree = reshape(adj_tree(:),N_alleles*N_alleles,[])';
    [~, ~, isomorphic] = unique(adj_tree, 'rows');
  
end
