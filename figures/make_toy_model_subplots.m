function make_toy_model_subplots(results_dir)

    plot_outdir = sprintf('%s/Figures/Simulation', results_dir);
    if (~exist(plot_outdir, 'dir'))
        mkdir(plot_outdir);
    end    
    plot_simulation_library_size();
    paper_print(sprintf('%s/PSingleton', plot_outdir));
    
    load(sprintf('%s/Simulation/ToyModel/Results.mat', results_dir));    
    plot_simulation_skew_and_depth(prolif_b, sample_fact, duplicate_prob, duplicate_prob_expansion);    
    paper_print(sprintf('%s/CollisionsObserved', plot_outdir));
    
    close all;

end