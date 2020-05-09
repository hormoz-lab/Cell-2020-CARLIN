function mt_mask = filter_out_genes_mt(genes)

    fprintf('Computing mitochondrial gene filter: ');
    
    mt_mask = uint16(find(startsWith(genes, 'mt-') | startsWith(genes, 'MT-')));
    
    fprintf('%d/%d\n', length(mt_mask), length(genes));
    
end