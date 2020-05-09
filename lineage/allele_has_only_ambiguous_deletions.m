function tf = allele_has_only_ambiguous_deletions(muts)
    
    if (isempty(muts))
        tf = false;
        return;
    end
    
    if (any(vertcat(muts.type)~='D'))
        tf = false;
        return;
    end

    N = size(muts,1);
    
    ref = CARLIN_def.getInstance;
    refseq = ref.seq.CARLIN;
    
    tf = true;
    for i = 1:size(muts,1)
        if (muts(i).loc_end < ref.width.CARLIN && refseq(muts(i).loc_start) ~= refseq(muts(i).loc_end+1) || ...
            muts(i).loc_start > 1 && refseq(muts(i).loc_end) ~= refseq(muts(i).loc_start-1))
            tf = false;
            return;
        end
    end
end