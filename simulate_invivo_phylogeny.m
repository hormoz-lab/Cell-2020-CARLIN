function simulate_invivo_phylogeny(analysis_dir, results_dir)    

    experiments = {{'InVivoPhylogeny/Rep1/Skin';    'InVivoPhylogeny/Rep2/Skin'};
                   {'InVivoPhylogeny/Rep1/Brain/L'; 'InVivoPhylogeny/Rep2/Brain/L'};
                   {'InVivoPhylogeny/Rep1/Brain/R'; 'InVivoPhylogeny/Rep2/Brain/R'};
                   {'InVivoPhylogeny/Rep2/Gr/L'};
                   {'InVivoPhylogeny/Rep2/Gr/R'};
                   {'InVivoPhylogeny/Rep2/BCell/L'};
                   {'InVivoPhylogeny/Rep2/BCell/R'};
                   {'InVivoPhylogeny/Rep1/MPP/L'};
                   {'InVivoPhylogeny/Rep1/MPP/R'};
                   {'InVivoPhylogeny/Rep1/HSC/L'};
                   {'InVivoPhylogeny/Rep1/HSC/R'};
                   {'InVivoPhylogeny/Rep1/Muscle/L';  'InVivoPhylogeny/Rep2/Muscle/L'};
                   {'InVivoPhylogeny/Rep1/Muscle/R';  'InVivoPhylogeny/Rep2/Muscle/R'};
                   {'InVivoPhylogeny/Rep1/Heart';     'InVivoPhylogeny/Rep2/Heart'};
                   {'InVivoPhylogeny/Rep1/Lung/L';    'InVivoPhylogeny/Rep2/Lung/L'};
                   {'InVivoPhylogeny/Rep1/Lung/R';    'InVivoPhylogeny/Rep2/Lung/R'};
                   {'InVivoPhylogeny/Rep1/Liver/1';   'InVivoPhylogeny/Rep1/Liver/2';   'InVivoPhylogeny/Rep2/Liver/1'; 'InVivoPhylogeny/Rep2/Liver/2'};                   
                   {'InVivoPhylogeny/Rep1/Intestine'; 'InVivoPhylogeny/Rep2/Intestine'};
                   {'InVivoPhylogeny/Rep1/Colon';     'InVivoPhylogeny/Rep2/Colon'}
                   };

    tissue_labels = {'Skin';
                     'BrainL';
                     'BrainR';
                     'GrL';
                     'GrR';
                     'BCellL';
                     'BCellR';
                     'MPPL';
                     'MPPR';
                     'HSCL';
                     'HSCR';
                     'MuscleL';
                     'MuscleR';
                     'Heart';
                     'LungL';
                     'LungR';
                     'Liver';
                     'Intestine';
                     'Colon'};

    for i = 1:length(experiments)
        for j = 1:length(experiments{i})
            temp{j} = load(sprintf('%s/%s/Summary.mat', analysis_dir, experiments{i}{j}), 'summary');
            temp{j} = temp{j}.summary;
        end
        samples{i} = ExperimentSummary.FromMerge(vertcat(temp{:}));
        clear temp;
    end
    samples = samples';
    filter_params = struct('epval', 0.05, 'cpval', 0, 'pos', false, 'neg', false, 'ambiguous', false, 'potential', 0);
    rng(1989);
    
    outdir = sprintf('%s/Trees/InVivoPhylogeny', results_dir);
    tissue_reconstruction(samples, outdir, results_dir, filter_params);
    save(sprintf('%s/Lineage.mat', outdir), 'tissue_labels', '-append');
end
