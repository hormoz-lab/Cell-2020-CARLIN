function [pval, tstat, df] = ttest_from_bootstrap(pos, neg, tail, N_pos, N_neg)

    if (nargin == 4)
        N_neg = N_pos;
    end

    var_pos = var(pos);
    var_neg = var(neg);
    
    mu_pos = mean(pos);
    mu_neg = mean(neg);
    
    df = ((var_pos+var_neg)^2)/(var_pos^2/N_pos+var_neg^2/N_neg);
    
    tstat = (mu_pos-mu_neg) / sqrt(var_pos+var_neg);
    
    if (strcmp(tail, 'right'))
        pval = 1-tcdf(tstat, df);
    else
        pval = tcdf(tstat, df);
    end
end