function [counts, barcodes, genes, platform] = load_count_matrix(input_count_file)

    if (endsWith(input_count_file, '.h5'))
        [counts, barcodes, genes] = Transcriptome.load_10x_h5_matrix(input_count_file);
        platform = '10X';
    elseif (endsWith(input_count_file, '.counts.tsv'))        
         [counts, barcodes, genes] = Transcriptome.load_indrops_matrix(input_count_file);
         platform = 'InDrops';
    else
        error('Unrecognized count matrix format\n');
    end
    
end