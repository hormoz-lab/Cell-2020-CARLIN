function [sample_sheet, samples_to_run] = get_experiment_metadata(experiment_series)

    [folder, ~, ~] = fileparts(mfilename('fullpath'));    
    sample_sheet = readtable(sprintf('%s/sample_sheet.csv', folder));    
    samples_to_run = find(ismember(sample_sheet.Series , experiment_series))';
    
end
