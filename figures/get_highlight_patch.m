function boundary_points = get_highlight_patch(x, y, candidate_set)    
    
    % Remove outliers    
    mu_x = median(x(candidate_set));
    mu_y = median(y(candidate_set));    
    diff_x = x(candidate_set)-mu_x;
    diff_y = y(candidate_set)-mu_y;
    diff_dist = diff_x.^2+diff_y.^2;
    prox_mask = find(diff_dist < percentile(diff_dist,95));
    prox_mask = prox_mask(boundary(x(candidate_set(prox_mask)), y(candidate_set(prox_mask)), 1));
    boundary_points = candidate_set(prox_mask);
    
    % Retrieve CCW ordering
    [~, reorder] = sort(atan2d(diff_y(prox_mask), diff_x(prox_mask)));
    boundary_points = boundary_points(reorder);
    
end