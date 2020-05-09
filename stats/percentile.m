function val = percentile(x, p)

    if (isempty(x))
        val = NaN;
        return;
    end
    if (length(x)==1)
        val = x;
        return;
    end
    
    assert(p >= 0 & p <= 100);
    
    x = sort(x);
    
    if (p == 0)
        val = x(1);
        return;
    end    
    if (p == 100)
        val = x(end);
        return;
    end
    
    ind = (length(x)-1)*(p/100.0)+1;
    frac = ind-fix(ind);
    ind = fix(ind);    
    val =  x(ind)+(x(ind+1)-x(ind))*frac;
    
end