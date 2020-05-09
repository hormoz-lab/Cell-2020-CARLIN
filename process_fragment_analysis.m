function [xvals, CDF, fraction_edited] = process_fragment_analysis(T, trim_amount, min_length, max_length, L, envelope)

    dat = [T.BeginPoint_BasePairs-trim_amount T.EndPoint_BasePairs-trim_amount T.Area_BasePairs];
    dat(dat(:,1) < min_length,:) = [];
    template_bins = find(dat(:,1)>=L-envelope(1) & dat(:,2)<=L+envelope(2));
    if (isempty(template_bins))
        template_bins = find(dat(:,1)>=L-envelope(1), 1, 'first');
    end
    if (isempty(template_bins))
        fraction_edited = 0;
    else
        fraction_edited = 1-sum(dat(template_bins,3))/sum(dat(:,3));
    end
    
    xvals = reshape(dat(:,1:2)', [numel(dat(:,1:2)), 1]);    
    CDF   = zeros(size(xvals));
    CDF(2:2:end) = cumsum(dat(:,3));
    xvals = [min_length; xvals];
    CDF = [0; CDF];
    CDF(2:2:end) = CDF(1:2:end-1);
    last_ind = find(xvals <= max_length, 1, 'last');
    xvals = [xvals(1:last_ind,:); max_length];
    CDF   = [CDF(1:last_ind); CDF(end)];
    CDF   = CDF/CDF(end);
    assert(length(xvals) == length(CDF));
    assert(issorted(xvals) && issorted(CDF));

end