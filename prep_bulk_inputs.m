function [input_file, cfg_type, output_path, num_cells] = prep_bulk_inputs(sample_sheet, s, processed_path, analysis_path)

    input_file = sprintf('%s/%s/PE.assembled.fastq.gz', processed_path, sample_sheet.ProcessedPath{s});
    cfg_type = sample_sheet.Type{s};
    output_path = sprintf('%s/%s', analysis_path, sample_sheet.ProcessedPath{s});
    num_cells = sample_sheet.CellsInSample{s};
    
    if (num_cells ~= '?')        
        num_cells = str2num(num_cells);
        if (strcmp(cfg_type, 'BulkRNA'))
            num_cells = num_cells*10;
        end
    end
    
end