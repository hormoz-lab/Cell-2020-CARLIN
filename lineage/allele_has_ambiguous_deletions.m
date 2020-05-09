function tf = allele_has_ambiguous_deletions(muts)
      
    N = size(muts,1);
    
    ref = CARLIN_def.getInstance;
    refseq = ref.seq.CARLIN;
    
    tf = false;
    
    for i = 1:size(muts,1)
        if (muts(i).type=='D')
            if (muts(i).loc_end < ref.width.CARLIN && refseq(muts(i).loc_start) == refseq(muts(i).loc_end+1) || ...
                muts(i).loc_start > 1 && refseq(muts(i).loc_end) == refseq(muts(i).loc_start-1))
                tf = true;
                return;
            end
        end
    end
end