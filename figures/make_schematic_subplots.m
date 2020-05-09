function make_schematic_subplots(results_dir)

    outdir = sprintf('%s/Figures/Schematics', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    plot_template_sequence('All');
    paper_print(sprintf('%s/TemplateFull', outdir), 'PNG');
    
    plot_template_sequence('Repeat');    
    paper_print(sprintf('%s/TemplateRepeat', outdir), 'PNG');
    
    plot_template_sequence('Prefix');    
    paper_print(sprintf('%s/TemplatePrefix', outdir), 'PNG');
    
    plot_template_sequence('Postfix');    
    paper_print(sprintf('%s/TemplatePostfix', outdir), 'PNG');
    
    close all;
    
end