function [cc_mask, cc_uncorr_mask, rho, pval, cc_signature] = filter_out_genes_cc_corr(normalized, scaled, genes, corr_cutoff)

    cc_genes = {'Ccnb1'; 'Plk1'; 'Cdc20'; 'Aurka'; 'Cenpf'; 'Cenpa'; 'Ccnb2'; 'Birc5'; 'Bub1'; 'Bub1b'; 'Ccna2'; 'Cks2'; 'E2f5'; 'Cdkn2d'};
    
    fprintf('%d genes marked as CC markers\n', length(cc_genes));
    
    [is, where] = ismember(cc_genes, genes);
    assert(all(is));
    cc_mask = uint16(where);
    
    cc_signature = nansum(scaled(:,cc_mask),2);
    
    [rho, pval] = corr(normalized, cc_signature);
    
    % This is just for subselecting genes for SPRING, so find stuff that
    % doesn't have high correlations, and let SPRING do whatever error
    % checking it normally does.
    cc_uncorr_mask = uint16(find(rho < corr_cutoff | isnan(rho)));
    
    fprintf('%d/%d genes correlated with cell cycle markers\n', length(genes)-length(cc_uncorr_mask), length(genes));
    
end
    
    
    
    