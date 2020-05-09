function [stable_adj_list, stable_allele_at_node, N_stable_alleles, stable_depthX_subtrees, N_omnipresent_depthX_subtrees, ...
          depthX_subtree, depthX_subtree_instances, depthX_subtree_pop, depthX_subtree_clade_pop, is_allele_at_depthX, depthX_alleles] = ...          
          compute_stable_tree(N_sim, N_alleles, N_samples, path_to_leaf, allele_at_node, allele_breakdown_by_sample)
      
    max_depth = 3;
    path_to_leaf = cellfun(@(x,y) x(y(:,1+[1:max_depth])), allele_at_node, path_to_leaf, 'un', false);
    
    % 1. For each simulation, determine the identity of nodes along a path length of
    %    up to 3 from the root to each leaf, and mark which alleles show up 
    %    at these depths.
    
    fprintf('Generating Gen1+2+3 Progenitors\n');

    depthX_progenitor_of_allele = cellfun(@squeeze, num2cell(shiftdim(cat(3, path_to_leaf{:}), 1),[2 3]), 'un', false);
   
    is_allele_at_depthX = repmat({false(N_sim, N_alleles)}, [max_depth,1]);
    
    % We pad paths up to max_depth, by repeating the leaf identities, so
    % remove these. Additionally, if an allele appears as both leaf and
    % internal node, record it at its earliest depth.
    
    for i = 1:N_sim        
        [r, ~, v] = find(accumarray(depthX_progenitor_of_allele{1}(i,:)',1));
        is_allele_at_depthX{1}(i, r(v>0)) = true;
        [r, ~, v] = find(accumarray(depthX_progenitor_of_allele{2}(i,:)',1));
        is_allele_at_depthX{2}(i, r(v>0)) = true;    
        is_allele_at_depthX{2}(i, is_allele_at_depthX{1}(i,:)) = false;
        [r, ~, v] = find(accumarray(depthX_progenitor_of_allele{3}(i,:)',1));
        is_allele_at_depthX{3}(i, r(v>0)) = true;    
        is_allele_at_depthX{3}(i, is_allele_at_depthX{1}(i,:)) = false;
        is_allele_at_depthX{3}(i, is_allele_at_depthX{2}(i,:)) = false;
    end

    % 2. Find all tree trunks of depth 1-3.
    %
    % ... this is a clunky way of doing it. You already have all the
    % progenitors so just do something like unique(...'row') conditioned on
    % depth.

    fprintf('Generating Progenitor Subtrees\n');

    depthX_alleles = arrayfun(@(g) arrayfun(@(i) find(is_allele_at_depthX{g}(i,:)), ...
                                                      [1:N_sim]', 'un', false), [1:max_depth]', 'un', false);
    [is, where] = arrayfun(@(g) arrayfun(@(i) ismember(depthX_progenitor_of_allele{g}(i,:), depthX_alleles{g}{i}), ...
                                              [1:N_sim]', 'un', false), [1:max_depth]', 'un', false);

    is_progenitor_depthX_subtree{1} = cellfun(@(x) find(x), is{1}, 'un', false);
    [depthX_subtree{1}, ~, which_depthX_subtree{1}] = cellfun(@(w1, i) unique(w1(i)'), ...
                                                                where{1}, is_progenitor_depthX_subtree{1}, 'un', false);
    depthX_subtree{1} = cellfun(@(w, a1) [a1(w(:,1))'], depthX_subtree{1}, ...
                                        depthX_alleles{1}, 'un', false);
    
    is_progenitor_depthX_subtree{2} = cellfun(@(x,y) find(x&y), is{1}, is{2}, 'un', false);
    [depthX_subtree{2}, ~, which_depthX_subtree{2}] = cellfun(@(w1, w2, i) unique([w1(i)' w2(i)'], 'rows'), ...
                                                                where{1}, where{2}, is_progenitor_depthX_subtree{2}, 'un', false);
    depthX_subtree{2} = cellfun(@(w, a1, a2) [a1(w(:,1))' a2(w(:,2))'], depthX_subtree{2}, ...
                                  depthX_alleles{1}, depthX_alleles{2}, 'un', false);
    
    is_progenitor_depthX_subtree{3} = cellfun(@(x,y,z) find(x&y&z), is{1}, is{2}, is{3}, 'un', false);
    [depthX_subtree{3}, ~, which_depthX_subtree{3}] = cellfun(@(w1, w2, w3, i) unique([w1(i)' w2(i)' w3(i)'], 'rows'), ...
                                                                where{1}, where{2}, where{3}, is_progenitor_depthX_subtree{3}, 'un', false);
    depthX_subtree{3} = cellfun(@(w, a1, a2, a3) [a1(w(:,1))' a2(w(:,2))' a3(w(:,3))'], depthX_subtree{3}, ...
                                  depthX_alleles{1}, depthX_alleles{2}, depthX_alleles{3}, 'un', false);

    N_premerge = arrayfun(@(g) cellfun(@(x) size(x,1), depthX_subtree{g}), [1:max_depth]', 'un', false);
    which_sim = arrayfun(@(g) repelem([1:N_sim]', N_premerge{g}), [1:max_depth], 'un', false);
    [depthX_subtree, ~, temp] = arrayfun(@(g) unique(vertcat(depthX_subtree{g}{:}), 'rows'), [1:max_depth]', 'un', false);    
    which_depthX_subtree = arrayfun(@(g) cellfun(@(x,y) x(y), mat2cell(temp{g}, N_premerge{g}), ...
                                                 which_depthX_subtree{g}, 'un', false), [1:max_depth]', 'un', false);
    
    % 3. For all subtrees of depth 1-3, find the distribution of 
    %    subclade populations.
    
    fprintf('Generating Gen1-3 Subclade Population Distributions\n');
    
    depthX_subtree_instances = arrayfun(@(g) accumarray(temp{g},1), [1:max_depth]', 'un', false);
    depthX_subtree_pop = arrayfun(@(g) arrayfun(@(v) zeros(v, N_samples), depthX_subtree_instances{g}, 'un', false), [1:max_depth]', 'un', false);    
    checksum_counter = cellfun(@(x) zeros(size(x)), depthX_subtree_instances, 'un', false);
    
    for g = 1:max_depth
        for t = 1:length(depthX_subtree{g})            
            for i = which_sim{g}(temp{g}==t)'                
                checksum_counter{g}(t) = checksum_counter{g}(t)+1;
                depthX_subtree_pop{g}{t}(checksum_counter{g}(t),:) = ...
                    sum(allele_breakdown_by_sample(is_progenitor_depthX_subtree{g}{i}(which_depthX_subtree{g}{i}==t),:),1);
            end
        end
    end

    assert(all(cellfun(@(x,y) isequal(x,y), depthX_subtree_instances, checksum_counter)));
    clear checksum_counter

    fprintf('Generating Progenitor Stats\n');

    depthX_subtree_clade_pop = arrayfun(@(g) struct('min', cellfun(@(x) min(sum(x,2)), depthX_subtree_pop{g}), ...
                                                    'max', cellfun(@(x) max(sum(x,2)), depthX_subtree_pop{g}), ...
                                                    'median', cellfun(@(x) median(sum(x,2)), depthX_subtree_pop{g})), [1:max_depth]', 'un', false);

    [~, reorder] = arrayfun(@(g) sortrows([depthX_subtree_instances{g}, depthX_subtree_clade_pop{g}.median], ...
                                          [1, 2], {'descend'; 'descend'}), [1:max_depth]', 'un', false);
    depthX_subtree = arrayfun(@(g) depthX_subtree{g}(reorder{g},:), [1:max_depth]', 'un', false);
    depthX_subtree_instances = arrayfun(@(g) depthX_subtree_instances{g}(reorder{g}), [1:max_depth]', 'un', false);
    depthX_subtree_pop       = arrayfun(@(g) depthX_subtree_pop{g}(reorder{g},:), [1:max_depth]', 'un', false);
    depthX_subtree_clade_pop = arrayfun(@(g) struct('min', depthX_subtree_clade_pop{g}.min(reorder{g}), ...
                                                    'max', depthX_subtree_clade_pop{g}.max(reorder{g}), ...
                                                    'median', depthX_subtree_clade_pop{g}.median(reorder{g})), [1:max_depth]', 'un', false);
    
    N_omnipresent_depthX_subtrees = arrayfun(@(g) find(depthX_subtree_instances{g}==N_sim, 1, 'last'), [1:max_depth]', 'un', false);
    stable_depthX_subtrees = arrayfun(@(g) depthX_subtree{g}(1:N_omnipresent_depthX_subtrees{g},:), [1:max_depth]', 'un', false);
    
    [stable_adj_list, stable_allele_at_node, N_stable_alleles] ...
            = generate_adj_list_from_subtrees(stable_depthX_subtrees, depthX_subtree_clade_pop, 'min', 1);   
   
end
