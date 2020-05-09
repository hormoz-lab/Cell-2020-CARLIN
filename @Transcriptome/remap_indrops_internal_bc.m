function remap_indrops_internal_bc(int_bc_file, abundant_bc_file, outfile)

    assert(exist(int_bc_file, 'file')==2, sprintf('Invalid path to InDrops internal barcode file: %s\n', int_bc_file));
    assert(exist(abundant_bc_file, 'file')==2, sprintf('Invalid path to InDrops abundant barcodes file: %s\n', abundant_bc_file));
    
    fid = fopen(int_bc_file);
    int_bc = textscan(fid, '%s');
    int_bc = int_bc{1};
    fclose(fid);
    
    fid = fopen(abundant_bc_file);
    abundant_bc = textscan(fid, '%s', 'Delimiter', ',');
    abundant_bc = abundant_bc{1};
    abundant_bc = reshape(abundant_bc, 3, [])';
    fclose(fid);
    
    abundant_bc(:,1) = strrep(abundant_bc(:,1), '-', '');
    
    [is, where] = ismember(int_bc, abundant_bc(:,2));
    assert(all(is));
    
    fid = fopen(outfile, 'wt');
    fprintf(fid, '%s\n', abundant_bc{where,1});
    fclose(fid);
    
end