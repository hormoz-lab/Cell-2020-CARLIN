function missing_gene = filter_out_genes_missing(A)

    fprintf('Computing missing gene filter: ');
    
    missing_gene = uint16(find(sum(A,1) == 0)');
    
    fprintf('%d/%d\n', length(missing_gene), size(A,2));
    
end