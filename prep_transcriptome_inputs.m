function [input_file, output_path] = prep_transcriptome_inputs(sample_sheet, s, processed_path, analysis_path)

    cfg_type = sample_sheet.Type{s};
    if (strcmp(cfg_type, 'CR2'))
        filename = 'filtered_gene_bc_matrices_h5.h5';
    elseif (strcmp(cfg_type, 'CR3'))
        filename = 'filtered_feature_bc_matrix.h5';
    end
    input_file = sprintf('%s/%s/%s', processed_path, sample_sheet.ProcessedPath{s}, filename);
    output_path = sprintf('%s/%s', analysis_path, sample_sheet.ProcessedPath{s});
    
end