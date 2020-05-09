function [z, pval] = two_proportion_ztest(np1, nq1, np2, nq2)
    
    N1 = np1+nq1;
    N2 = np2+nq2;

    p1 = np1/N1;
    p2 = np2/N2;
    p = (np1+np2)/(N1+N2);
    
    z = (p1-p2)/sqrt(p*(1-p)*(1/N1+1/N2));
    pval = normcdf(z, 'upper');

end