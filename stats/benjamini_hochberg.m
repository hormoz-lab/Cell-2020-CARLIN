function p_crit = benjamini_hochberg(pval, fdr_level)

    pval = sort(pval);
    n = length(pval);    
    qvalues = [1:n]'/n*fdr_level;
    crit_ind = find(pval <= qvalues, 1, 'last');
    if (isempty(crit_ind))
        p_crit = 0;
    else
        p_crit = pval(crit_ind);
    end
    
end