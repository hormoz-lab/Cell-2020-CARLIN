function plot_mutations_after_induction_round(yvals, series_names)

    N_series = length(series_names);
    assert(size(yvals,2)==N_series);
    
    fig_height = 4;
    fig_width  = 4.275;    

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
    
    plot(yvals, 'LineWidth', 1.0); 
    axis tight;
    ylim([0 1]);
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(gca, 'xtick', [1:size(yvals,1)], 'xticklabel', cellstr(num2str([1:size(yvals,1)]')));
    xlabel('Number of Mutations', 'FontSize', 6);
    
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    ylabel('Fraction Alleles in Induction Round', 'FontSize', 6);
    
    set(gca, 'LineWidth', 1.0);    
    set(gca, 'TickDir', 'out');
    
    box off;
    
    legend(series_names, 'Location', 'Northeast', 'FontSize', 5);
    legend('boxoff');

end