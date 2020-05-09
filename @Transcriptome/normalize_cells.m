function A = normalize_cells(A, exclude_frac)

    fprintf('Normalizing cell expression levels...\n');

    if (nargin < 2)
        exclude_frac = 1.0;
    end
    
    domineering_genes = sum(spdiags(1./sum(A,2), 0, size(A,1), size(A,1))*A > exclude_frac) > 0;    
    tot_UMI = sum(A(:,~domineering_genes),2);
    A = spdiags(mean(tot_UMI)./tot_UMI, 0, size(A,1), size(A,1))*A;
    
end
    
   




