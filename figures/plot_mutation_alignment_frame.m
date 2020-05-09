function sp = plot_mutation_alignment_frame(dat, title_str, series_names, show_legend)
    
    if (nargin < 4)
        show_legend = false;
    end
    
    ref = CARLIN_def.getInstance;
    W = ref.width.segment+ref.width.pam;
    
    assert(size(dat,1) == W && size(dat,2)==4);
    
    fig_width  = 5.8;
    fig_height = 3.0;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);       
    sp = subplot(1,1,1);
    
    max_y = max(dat(:));
    
    hold on;
    
    patch([0.5 0.5 ref.width.consite+0.5 ref.width.consite+0.5], [0 max_y max_y 0], [1 1 1]*ref.alpha.consite, 'FaceAlpha', ref.alpha.overlay);
    patch(ref.width.consite+[0.5 0.5 ref.width.cutsite+0.5 ref.width.cutsite+0.5], [0 max_y max_y 0], [1 1 1]*ref.alpha.cutsite, 'FaceAlpha', ref.alpha.overlay);
    patch(ref.width.segment+[0.5 0.5 ref.width.pam+0.5 ref.width.pam+0.5], [0 max_y max_y 0], [1 1 1]*ref.alpha.pam, 'FaceAlpha', ref.alpha.overlay);        
    bar([1:W]-0.375, dat(:,1), 'b', 'BarWidth', 0.25, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    bar([1:W]-0.125, dat(:,2), 'b', 'BarWidth', 0.25, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    bar([1:W]+0.125, dat(:,3), 'r', 'BarWidth', 0.25, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    bar([1:W]+0.375, dat(:,4), 'r', 'BarWidth', 0.25, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
   
    hold off;
    
    axis tight;
    
    max_expo = floor(log10(max_y));
    set(get(gca, 'YAxis'), 'Exponent', max_expo);
    
    xticks([1:W]);
    xticklabels(cellstr(num2str([1:ref.width.consite 1:ref.width.cutsite 1:ref.width.pam]')));
    xtickangle(90);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('       Conserved Site              Cutsite        PAM+Linker', 'FontSize', 6);
    ylabel('Frequency In Alleles', 'FontSize', 6);    
        
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
        leg_width  = 2.0;
        lgd.Position = [left_margin+tight_margin fig_height-top_margin-leg_height-4*tight_margin leg_width leg_height];
        legend('boxoff');
    end
    
end