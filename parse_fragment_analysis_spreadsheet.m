function T = parse_fragment_analysis_spreadsheet(file, sheet)
    assert(exist(file, 'file')==2);
    T = readtable(file, 'Sheet', sheet);    
    T.Properties.VariableNames = cellfun(@(x) strip(x, 'right', '_'), T.Properties.VariableNames, 'un', false);
end
