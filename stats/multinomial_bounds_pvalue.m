% Corrado, C.J. 2011. "The exact distribution of the maximum, minimum and the range of 
% Multinomial/Dirichlet and Multivariate Hypergeometric frequencies." Statistics and computing 21 (3): 349-359.

function [pval, out, Q, count_mask] = multinomial_bounds_pvalue(N, p, min_count, max_count)
    
    assert((sum(p)-1)<1e-10);
    M = length(p);
    
    Q = cell(M,1);
    
    binom_coeff = zeros(N+1,N+1);
    binom_coeff(:,1) = 1;
    
    for n = 2:N+1
        for k = 2:n
            binom_coeff(n,k) = binom_coeff(n-1,k)+binom_coeff(n-1,k-1);
        end
    end
    
    binom_coeff = fliplr(flipud(binom_coeff));
    exp_coeff = arrayfun(@(k) diag(k*ones(N-k+1, 1), k), [1:N]', 'un', false);
    exp_coeff = sum(cat(3,reshape(horzcat(exp_coeff{:}), [N+1, N+1, N])),3);    
    pstar = p./cumsum(p, 'rev');
    
    Q(1:M-1) = arrayfun(@(k) binom_coeff.*((pstar(k)).^exp_coeff).*((1-pstar(k)).^([N:-1:0])), [1:M-1]', 'un', false);           
    Q{M} = zeros(N+1,N+1);
    Q{M}(:,end) = 1;
    
    count_mask = arrayfun(@(k) triu(true(N+1), min_count(k)).*tril(true(N+1), max_count(k)), [1:M]', 'un', false);
    Q = arrayfun(@(k) Q{k}.*count_mask{k}, [1:M]', 'un', false);
    
    out = eye(N+1);
    for k = 1:M
        out = out*Q{k};
    end
    
    pval = out(1,end);
    
end