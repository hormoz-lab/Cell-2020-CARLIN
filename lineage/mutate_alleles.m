function [relationship, mutations, type, preserved] = mutate_alleles(parent_seq, child_seq)
    
    assert(isa(parent_seq, 'AlignedSEQ') && isa(child_seq, 'AlignedSEQ'));
    
    % Get mutation list for each aligned sequence
    parent_mut = Mutation.identify_Cas9_events(parent_seq);
    child_mut  = Mutation.identify_Cas9_events( child_seq);
    
    [relationship, mutations, type, preserved] = mutate_mut_list(parent_mut, child_mut);    
end