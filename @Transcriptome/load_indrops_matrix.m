function [counts, barcodes, genes] = load_indrops_matrix(count_file)

    assert(exist(count_file, 'file')==2, sprintf('Invalid path to InDrops counts file: %s\n', count_file));
    
    fprintf('Loading InDrops counts matrix: %s\n', count_file);

    fid = fopen(count_file);
    genes = fgets(fid);
    genes = strsplit(genes)';
    genes = genes(2:end-1);    
    barcodes = textscan(fid, '%s %*[^\n]');
    barcodes = barcodes{1};
    fclose(fid);

    counts = sparse(dlmread(count_file, '', 1, 1));
    
    assert(size(counts,1) == size(barcodes,1), 'Number of barcodes and rows in count matrix do not match');
    assert(size(counts,2) == size(genes,1), 'Number of genes and columns in count matrix do not match');
    
end