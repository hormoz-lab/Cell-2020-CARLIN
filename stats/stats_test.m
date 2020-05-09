function tests = stats_test
    tests = functiontests(localfunctions);
end


% Corrado, C.J. 2011. "The exact distribution of the maximum, minimum and the range of 
% Multinomial/Dirichlet and Multivariate Hypergeometric frequencies." Statistics and computing 21 (3): 349-359.

function test_bounds_probability(testCase)

    % Example from Section 2 of Corrado

    [p, ~, Q, ~] = multinomial_bounds_pvalue(3, ones(3,1)/3, [0 0 0], [3 3 3]);
    verifyLessThanOrEqual(testCase, abs(p-1), 1e-6);
    verifyLessThanOrEqual(testCase, max(abs(Q{1}(1,:)-[0.296296 0.444444 0.222222 0.037037])), 1e-6);
    verifyLessThanOrEqual(testCase, max(max(abs(Q{2}-[0.125 0.375 0.375 0.125;
                                                      0     0.25  0.5   0.25; 
                                                      0     0     0.5   0.5; 
                                                      0     0     0     1]))), 1e-6);
    verifyLessThanOrEqual(testCase, max(abs(Q{3}(:,end)-ones(4,1))), 1e-6);

    % Example from Section 3 of Corrado (max, c=2)

    [p, ~, Q, ~] = multinomial_bounds_pvalue(3, ones(3,1)/3, [0 0 0], [2 2 2]);
    verifyLessThanOrEqual(testCase, abs(p-0.888889), 1e-6);
    verifyLessThanOrEqual(testCase, max(abs(Q{1}(1,:)-[0.296296 0.444444 0.222222 0])), 1e-6);
    verifyLessThanOrEqual(testCase, max(max(abs(Q{2}-[0.125 0.375 0.375 0;
                                                      0     0.25  0.5   0.25; 
                                                      0     0     0.5   0.5; 
                                                      0     0     0     1]))), 1e-6);
    verifyLessThanOrEqual(testCase, max(abs(Q{3}(:,end)-[0; 1; 1; 1])), 1e-6);

    % Example from Section 3 of Corrado (min, c=1)

    [p, ~, Q, ~] = multinomial_bounds_pvalue(3, ones(3,1)/3, [1 1 1], [3 3 3]);
    verifyLessThanOrEqual(testCase, abs(p-0.222222), 1e-6);
    verifyLessThanOrEqual(testCase, max(abs(Q{1}(1,:)-[0 0.444444 0.222222 0.037037])), 1e-6);
    verifyLessThanOrEqual(testCase, max(max(abs(Q{2}-[0    0.375 0.375 0.125;
                                                      0    0     0.5   0.25; 
                                                      0    0     0     0.5; 
                                                      0    0     0     0]))), 1e-6);
    verifyLessThanOrEqual(testCase, max(abs(Q{3}(:,end)-[1; 1; 1; 0])), 1e-6);

    % Example from Section 3 of Corrado (c <= 19)

    p = multinomial_bounds_pvalue(500, ones(50,1)/50, zeros(50,1), 19*ones(50,1));
    verifyLessThanOrEqual(testCase, abs(p-0.852727), 1e-6);

    % Example from Section 3 of Corrado (c >= 3) - I think there's a typo in the paper, should read less than 3, not 4

    p = multinomial_bounds_pvalue(500, ones(50,1)/50, 3*ones(50,1), 500*ones(50,1));
    verifyLessThanOrEqual(testCase, abs(p-0.877373), 1e-6);

    % Example from Section 3 of Corrado (3 <= c <= 19) - I think there's a typo in the paper, should read less than 3, not 4
    p = multinomial_bounds_pvalue(500, ones(50,1)/50, 3*ones(50,1), 19*ones(50,1));
    verifyLessThanOrEqual(testCase, abs(p-0.750895), 1e-6);

end

function test_range_pvalue(testCase)
    
    % From table 1 of Corrado
    
    p = arrayfun(@(r) multinomial_range_pvalue(25, ones(4,1)/4, r), [3:10]');
    act = [0.1046 0.2450 0.4184 0.5966 0.7483 0.8599 0.9296 0.9681]';
    verifyLessThanOrEqual(testCase, max(abs(p-act)), 2e-4);
    
    p = arrayfun(@(r) multinomial_range_pvalue(50, ones(8,1)/8, r), [3:10]');
    act = [0.0037 0.0256 0.0933 0.2279 0.4149 0.6114 0.7743 0.8845]';
    verifyLessThanOrEqual(testCase, max(abs(p-act)), 1e-4);
    
    p = arrayfun(@(r) multinomial_range_pvalue(125, ones(25,1)/25, r), [3:10]');
    act = [0.0000 0.0000 0.0009 0.0103 0.0736 0.2440 0.4994 0.7348]';
    verifyLessThanOrEqual(testCase, max(abs(p-act)), 1e-4);

end



