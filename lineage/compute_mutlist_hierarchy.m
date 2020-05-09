function [mut_matrix, relationship, mutations, type, preserved] = ...
          compute_mutlist_hierarchy(parent_muts, child_muts)
      
    N1 = size(parent_muts,1);
    
    if (nargin == 2)
        N2 = size(child_muts,1);        
    else
        N2 = N1;
        child_muts = parent_must;
    end
    
    relationship = cell(N1,N2);
    mutations    = cell(N1,N2);
    type         = cell(N1,N2);
    preserved    = cell(N1,N2);
    
    mut_matrix = inf(N1,N2);
    
    parfor i = 1:N1
        i
        for j = 1:N2
            [i j]
            [relationship{i,j}, mutations{i,j}, type{i,j}, preserved{i,j}] = mutate_mut_list(parent_muts{i}, child_muts{j});
            if ~(strcmp(relationship{i,j}, 'Incompatible'))
                mut_matrix(i,j) = length(mutations{i,j});
            end            
        end
    end
end