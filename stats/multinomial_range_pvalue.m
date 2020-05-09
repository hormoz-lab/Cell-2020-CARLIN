% Corrado, C.J. 2011. "The exact distribution of the maximum, minimum and the range of 
% Multinomial/Dirichlet and Multivariate Hypergeometric frequencies." Statistics and computing 21 (3): 349-359.

function [pval, out, Q, count_mask] = multinomial_range_pvalue(N, p, r)

    if (r > N+1)
        r = N+1;
    end
    
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
    
    t1 = repmat({eye(N+1)}, [N-r+2, 1]);    
    for h = 0:(N-r+1)
        count_mask = triu(true(N+1), h).*tril(true(N+1), h+r-1);
        Qh = arrayfun(@(k) Q{k}.*count_mask, [1:M]', 'un', false);
        for k = 1:M
            t1{h+1} = t1{h+1}*Qh{k};
        end
    end
    
    if (N-r+1 > 0)
        t2 = repmat({eye(N+1)}, [N-r+1, 1]);
        for h = 0:(N-r)
            count_mask = triu(true(N+1), h+1).*tril(true(N+1), h+r-1);
            Qh = arrayfun(@(k) Q{k}.*count_mask, [1:M]', 'un', false);
            for k = 1:M
                t2{h+1} = t2{h+1}*Qh{k};
            end
        end        
        out = sum(cat(3, t1{:}),3)-sum(cat(3, t2{:}),3);    
    else
        out = sum(cat(3, t1{:}),2);
    end
    pval = out(1,end);
    
end