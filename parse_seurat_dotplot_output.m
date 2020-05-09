function out = parse_seurat_dotplot_output(read_path, marker_genes, louvain_cluster)

    out = readtable(sprintf('%s/DotPlot.csv', read_path));
    out = out(:,[2:end]);
    out.Properties.VariableNames = {'AvgExp'; 'PctExp'; 'Gene'; 'Louvain'; 'AvgExpScaled'};
    out.Louvain = cellfun(@str2num, out.Louvain);
    [is, out.WhichMarkerGenes] = ismember(out.Gene, marker_genes);
    assert(all(is));
    [is, out.WhichLouvainCluster] = ismember(out.Louvain, louvain_cluster);
    assert(all(is));
    
end
        