function [input_file, cfg_type, output_path, CB_file] = prep_SC_inputs(sample_sheet, s, processed_path, analysis_path)

    cfg_type = sample_sheet.Type{s};
    input_path = sprintf('%s/%s', processed_path, regexprep(sample_sheet.ProcessedPath{s}, 'Amplicon.*', 'Amplicon'));
    input_file = cellfun(@(x) {sprintf('%s/%s_R1_001.fastq.gz', input_path, x), ...
                               sprintf('%s/%s_R2_001.fastq.gz', input_path, x)}, strsplit(sample_sheet.RawFile{s}, ';')', 'un', false);
    input_file = vertcat(input_file{:});
    output_path = sprintf('%s/%s', analysis_path, sample_sheet.ProcessedPath{s});
    CB_file = sprintf('%s/filtered_barcodes_umi_mt.txt', regexprep(output_path, 'Amplicon.*', 'Transcriptome'));
    
end