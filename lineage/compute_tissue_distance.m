function tissue_distance = compute_tissue_distance(adj_list, allele_breakdown_by_sample)

    [N_alleles, N_samples] = size(allele_breakdown_by_sample);
    N_nodes = length(adj_list);
    
    tissue_distance_num = zeros(N_samples, N_samples);
    tissue_distance_den = zeros(N_samples, N_samples);
    
    allele_breakdown_by_sample = padarray(allele_breakdown_by_sample, [N_nodes-N_alleles, 0], 0, 'post');    
    
    for k = 2:N_alleles            
        cur = k;
        while (adj_list(cur) ~= 0)
            cur = adj_list(cur);
            allele_breakdown_by_sample(cur,:) = allele_breakdown_by_sample(cur,:) + allele_breakdown_by_sample(k,:);            
        end
    end
    
    assert(sum(allele_breakdown_by_sample(N_alleles+1,:))==sum(sum(allele_breakdown_by_sample(2:N_alleles,:))));
    
    for i = 1:N_samples
        alleles_expressing_tissue = find(allele_breakdown_by_sample(2:N_alleles,i))+1;
        for a = alleles_expressing_tissue'
            tissues_yet_to_encounter = [1:N_samples];
            d = 0;
            cur = a;
            while(adj_list(cur) ~= 0 && ~isempty(tissues_yet_to_encounter))
                encountered_tissues = find(allele_breakdown_by_sample(a,:));
                if (~isequal(encountered_tissues,i))
                    tissues_to_update = intersect(tissues_yet_to_encounter, encountered_tissues);
                    for j = tissues_to_update
                        tissue_distance_num(i,j) = tissue_distance_num(i,j) + d*allele_breakdown_by_sample(a,i)*allele_breakdown_by_sample(cur,j);
                        tissue_distance_den(i,j) = tissue_distance_den(i,j) + allele_breakdown_by_sample(a,i)*allele_breakdown_by_sample(cur,j);
                    end                
                    d = d+1;
                    tissues_yet_to_encounter = setdiff(tissues_yet_to_encounter, tissues_to_update);
                end
                cur = adj_list(cur);
            end
            if (~isempty(tissues_yet_to_encounter))
                assert(all(allele_breakdown_by_sample(cur,tissues_yet_to_encounter)>0));
                for j = tissues_yet_to_encounter
                    tissue_distance_num(i,j) = tissue_distance_num(i,j) + d*allele_breakdown_by_sample(a,i)*allele_breakdown_by_sample(cur,j);
                    tissue_distance_den(i,j) = tissue_distance_den(i,j) + allele_breakdown_by_sample(a,i)*allele_breakdown_by_sample(cur,j);
                end
            end 
        end
    end
    
    tissue_distance = tissue_distance_num./tissue_distance_den;
    tissue_distance(isnan(tissue_distance)) = inf;
    
    tissue_distance = tissue_distance+tissue_distance';
    
end     
