function A = normalize_genes(A, recentre)

    fprintf('Computing gene z-scores...\n');

    if (nargin < 2)
        recentre = true;
    end
    
    if (recentre)
       A = full((A-mean(A,1))./std(A,0,1));
    else
       A = A * spdiags(1./std(A,0,1)', 0, size(A,2), size(A,2));
    end
    
end
    
   