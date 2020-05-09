function [mut_matrix, relationship, mutations, type, preserved, parent_muts, child_muts] = ...
          compute_allele_hierarchy(aligned_seqs1, aligned_seqs2)
  
    N1 = size(aligned_seqs1,1);    
    parent_muts = cell(N1,1);
    parfor i = 1:N1
        parent_muts{i} = Mutation.identify_Cas9_events(aligned_seqs1{i});
    end
    
    if (nargin == 2)
        N2 = size(aligned_seqs2,1);
        child_muts  = cell(N2,1);
        parfor i = 1:N2
            child_muts{i} = Mutation.identify_Cas9_events(aligned_seqs2{i});
        end
    else
        N2 = N1;
        child_muts = parent_muts;
    end
    
    [mut_matrix, relationship, mutations, type, preserved] = compute_mutlist_hierarchy(parent_muts, child_muts);
 
end