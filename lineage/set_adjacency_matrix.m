function adj = set_adjacency_matrix(mut_matrix, relationship, mutations, type, preserved, allele_freq)

    N_alleles = size(mut_matrix,1);
    self_map = logical(eye(N_alleles));
    
    adj = strcmp(relationship, 'Offspring');
    
    adj = adj & allele_freq > allele_freq';
    
    num_preserved = cellfun(@(x) sum(x), preserved);
    num_preserved(self_map) = 0;
    num_preserved(~adj) = 0;
    max_preserved = max(num_preserved,[],1);
    adj = adj & num_preserved==max_preserved & max_preserved > 0;
        
    mut_matrix(self_map) = inf;
    mut_matrix(~adj) = inf;
    adj = adj & mut_matrix==min(mut_matrix,[],1);
    
    adj(1,~any(adj,1)) = true;
    adj(self_map) = true;
    
end