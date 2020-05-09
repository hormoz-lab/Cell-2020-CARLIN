function [pval, N_min] = fate_bias_pvalue(summary, allele_breakdown_by_sample, sig_level)

    template_mask = cellfun(@(x) isequal(degap(x.get_seq), CARLIN_def.getInstance.seq.CARLIN), summary.alleles);
    
    allele_breakdown_by_sample = allele_breakdown_by_sample(~template_mask, :);
    
    N = sum(allele_breakdown_by_sample,2);
    
    % We should calculate probabilities omitting the row to test to be
    % strictly correct, but for some samples there are a few rows that
    % dominant. Leaving the row in for computing probabilities will make
    % the p-value more conservative, since the null distribution will be
    % weighted to look more like the sample we test against.
    
    p = sum(allele_breakdown_by_sample,1)/sum(N);
    assert(abs(sum(p)-1)<= 1e-10);
    
    pval = zeros(size(summary.alleles));    
    pval(~template_mask) = arrayfun(@(i) multinomial_range_pvalue(N(i), p, max(allele_breakdown_by_sample(i,:))-min(allele_breakdown_by_sample(i,:))), [1:size(N,1)]');
    pval = 1-pval;
    
    if (nargout > 1)
        N_min = 1;        
        while (multinomial_range_pvalue(N_min, p, N_min) < 1-sig_level)
            N_min = N_min+1;
        end
    end
    
end