function plot_simulation_library_size()
    
    fig_width  = 7.0;
    fig_height = 3.9;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    N_alleles = [150; 500; 1000; 5000; 10000; 20000; 44000];
    N_cells = round(logspace(2, 5, 10)');
    p_singleton = arrayfun(@(p) arrayfun(@(N) nchoosek(N,1)*p*(1-p)^(N-1)/...
                                              (1-nchoosek(N,0)*p^0*(1-p)^N), N_cells), 1./N_alleles, 'un', false);
    p_singleton = horzcat(p_singleton{:});            
    
    p = plot(N_cells, p_singleton, 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 1.0);
    set(p, {'MarkerFaceColor'}, get(p,'Color'));
    set(gca, 'xscale', 'log');
    set(get(gca, 'xaxis'), 'FontSize', 5);
    set(get(gca, 'yaxis'), 'FontSize', 5);    
    xlabel('Observed Cells', 'FontSize', 6);
    ylabel('Lower bound on p_{singleton}', 'FontSize', 6);
    axis tight;
    ylim([0 1]);
    title(sprintf('Uniquely Marking Allele Fraction'), ...
        'FontWeight', 'normal', 'FontSize', 6);
    lgd = legend(cellstr(num2str(N_alleles)), 'FontSize', 5, 'Location', 'EastOutside');
    title(lgd, 'Number of Alleles', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    box off;
    
    left_margin = 0.8;
    top_margin = 0.3;    
    bottom_margin = 0.8;
    right_margin = 0.2;
    leg_width = 2.0;    
    leg_height = 1.0;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin-leg_width fig_height-top_margin-bottom_margin]);

    lgd.Units = 'centimeters'; 
    lgd.Position = [fig_width-leg_width, fig_height/2, leg_width, leg_height];
    
end