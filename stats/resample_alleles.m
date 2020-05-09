function [resampled, clone_size_of_allele, clone_size_freqs, transcripts_per_clone_size] = resample_alleles(datavec, N_sample, N_trials)

    % Datavec are allele labels {j}>0 of samples {i} such that X_i = j,    
    
    % Draw N_sample times with replacement from datavec, and repeat over
    % N_trial bootstrap trials
    
    % resampled                  - [1 x N_trial] cell-array, each with [N_sample x 1] draws from datavec
    % clone_size_of_alleles      - population of each allele over N_trials
    % clone_size_freqs           - the number of clones of a given size > 0 over N_trials
    % transcripts_per_clone_size - the number of transcripts belonging to a clone of a given size > 0, over N_trials
    
    assert(all(datavec>0));
    
    % Maximum index of allele label
    N_alleles = max(datavec);
    
    % Number of unique clones in original sample
    N_clones = length(unique(datavec));

    resampled  = arrayfun(@(i) datasample(datavec, N_sample), [1:N_trials], 'un', false);
    
    if (nargout > 1)
        
        clone_size_of_allele = cellfun(@(x) accumarray(x, 1, [N_alleles, 1]), resampled, 'un', false);
        
        if (nargout > 2)
    
            clone_size_freqs = cellfun(@(x) accumarray(nonzeros(x), 1, [N_sample, 1]), clone_size_of_allele, 'un', false);            
            % Can't ever see more clones than existed in the original data            
            assert(all(cellfun(@(x) sum(x) <= N_clones, clone_size_freqs)));
            
            if (nargout > 3)
    
                transcripts_per_clone_size = cellfun(@(x) x.*[1:N_sample]', clone_size_freqs, 'un', false);                
                % All transcripts have to belong to clones of a given size                
                assert(all(cellfun(@(x) sum(x)==N_sample, transcripts_per_clone_size)));
                
            end
        end
    end
      
end