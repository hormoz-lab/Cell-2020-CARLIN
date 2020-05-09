% Orlitsky, Alon, Ananda Theertha Suresh, and Yihong Wu. 2016. "Optimal prediction of the number of unseen species." 
% Proceedings of the National Academy of Sciences 113 (47): 13283-13288.

function U = Orlitzky(m, freq_counts)

    [r, ~, v] = find(freq_counts);
    n = sum(r.*v);
    assert(m <= n*log(n));
    t = m/n;
    
    if (t < 1)
        U = -sum(((-t).^r).*v);
    else
        w = binopdf([1:r(end)]',floor(log(n*t^2/(t-1))/log(3)/2),2/(t+2));
        w = flipud(cumsum(flipud(w)));
        U = v.*(w(r).*((-t).^r));
        U(isnan(U)) = 0;
        U = -sum(U);
    end
end
