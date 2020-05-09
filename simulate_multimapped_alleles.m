function simulate_multimapped_alleles(results_dir)

    % Just pick a reasonable duplication probability
    % Since we are normalizing by this our choice doesn't really matter
    duplicate_prob = 0.05;
    N_alleles = round(linspace(150, 44000, 10)');
    syms N_cells
    eqn = arrayfun(@(p) (1-nchoosek(N_cells,0)*p^0*(1-p)^N_cells-nchoosek(N_cells,1)*p*(1-p)^(N_cells-1))/...
                        (1-nchoosek(N_cells,0)*p^0*(1-p)^N_cells)==duplicate_prob, 1./N_alleles, 'un', false);
    
    N_cells = cellfun(@(x) round(double(vpasolve(x, N_cells))), eqn);

    N_sim = 1000;
    prolif_b = [0; 1; 2; 5; 10; 20; 100];
    sample_fact = [0.5; 1.0; 5.0; 10; 100];

    % For a fixed dup_prob, the number of alleles doesn't matter (modulo edge
    % effects), so just pick one as an exemplar. If I pick a different
    % value, the sampling depth plot doesn't change.
    exemplar = 2; % N_alleles(exemplar) ~= 5e3

    duplicate_prob_expansion = zeros(length(prolif_b), length(sample_fact));

    allele_labels = arrayfun(@(k) datasample([1:N_alleles(exemplar)]', N_cells(exemplar), 'Replace', true), [1:N_sim]', 'un', false);
    label_multiplicity = cellfun(@(x) accumarray(x, 1, [N_alleles(exemplar), 1]), allele_labels, 'un', false);        
    fprintf('Target Prob: %f, Achieved Prog: %f\n', duplicate_prob, mean(cellfun(@(x) sum(x>1)/sum(x>0), label_multiplicity)));
    for j = 1:length(prolif_b)
        p = (1-[0:N_cells(exemplar)-1]'./(N_cells(exemplar)-1)).^prolif_b(j);
        p = cumsum(p)/sum(p);
        for k = 1:length(sample_fact)
            u = rand(round(sample_fact(k)*N_cells(exemplar)),N_sim);
            prog = quickinv(u,p);
            prog = cellfun(@(x,y) [x y(x)], num2cell(prog,1)', allele_labels, 'un', false);        
            label_collision = cellfun(@(x) unique(x, 'rows'), prog, 'un', false);
            label_collision = cellfun(@(x) accumarray(x(:,2), 1), label_collision, 'un', false);
            duplicate_prob_expansion(j,k) = mean(cellfun(@(x) sum(x>1)/sum(x>0), label_collision));
        end
    end                    
    
    outdir = sprintf('%s/Simulation/ToyModel', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    save(sprintf('%s/Results.mat', outdir));
    
end

function prog = quickinv(u, p)

    [N_draw, N_sim] = size(u);

    [u, order] = sort(u(:));
    c = 1;
    cur = zeros(length(p),1);
    for k = 1:length(p)
        while(u(c) <= p(k) && c < length(u))
            c = c+1;
        end
        cur(k) = c-1;
    end
    cur = [0; cur];
    prog = arrayfun(@(i) i*ones(cur(i+1)-cur(i),1), [1:length(cur)-1]', 'un', false);
    prog = vertcat(prog{:});    
    if (length(prog) == length(u)-1)
        prog = [prog; length(p)];
    end
    assert(length(u) == length(prog));
    
    [~, reorder] = sort(order);
    prog = prog(reorder);
    
    prog = reshape(prog, [N_draw, N_sim]);
end
