function sp = plot_mutation_alignment_count(dat, title_str, series_names, show_legend)

    if (nargin < 4)
        show_legend = false;
    end
    
    N_trunc = size(dat,1)-1;
    assert(size(dat,2)==4);
    
    fig_width  = 5.8;
    fig_height = 3.0;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);       
    sp = subplot(1,1,1);
    
    hold on;
    bar([0:N_trunc]-0.3, dat(:,1), 'b', 'BarWidth', 0.2, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    bar([0:N_trunc]-0.1, dat(:,2), 'b', 'BarWidth', 0.2, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    bar([0:N_trunc]+0.1, dat(:,3), 'r', 'BarWidth', 0.2, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    bar([0:N_trunc]+0.3, dat(:,4), 'r', 'BarWidth', 0.2, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    hold off;
    axis tight;
    
    max_y = ylim;
    max_y = max_y(2);
    max_expo = floor(log10(max_y));
    set(get(gca, 'YAxis'), 'Exponent', max_expo);
    
    xticks([0:N_trunc]);
    xticklabels([cellstr(num2str([0:N_trunc-1]')); {sprintf('>=%d', N_trunc)}]);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('Number of Mutations in Allele', 'FontSize', 6);
    ylabel('Number of Alleles', 'FontSize', 6);    
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    box off;
    
    left_margin = 0.7;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    tight_margin = 0.1;        
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-tight_margin fig_height-top_margin-bottom_margin]);
    
    if (show_legend)
        lgd = legend(series_names, 'FontSize', 5, 'Location', 'Northwest');
        lgd.Units = 'centimeters';
        leg_height = 1.0;
        leg_width  = 3.0;
        lgd.Position = [fig_width-leg_width-tight_margin fig_height-top_margin-leg_height-4*tight_margin leg_width leg_height];
        legend('boxoff');
    end
    
end