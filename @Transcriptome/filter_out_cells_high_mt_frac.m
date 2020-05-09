function [mt_mask, mt_frac_cutoff] = filter_out_cells_high_mt_frac(A, mt_genes, mt_frac_cutoff)

    fprintf('Filtering high mitochrondrial fraction cells: ');

    if (nargin < 3 || mt_frac_cutoff < 0)
        mt_frac_cutoff = 0.15;
    end
    
    mt_mask = full(sum(A(:,mt_genes),2)./sum(A,2) <= mt_frac_cutoff);   
    fprintf('%d/%d survive at MT frac cutoff <= %3.2g\n', sum(mt_mask), length(mt_mask), mt_frac_cutoff);
    
end
    
    
    
    