function sp = plot_phenotype_on_transcriptome(x, y, ceres, ceres_names, plot_order_bg_to_fg, title_str, body, double_column, legend_reorder)
        
    if (body)
        fig_width = 4.6;
        fig_height = 3.5;
        double_column = false;
    else
        fig_width  = 8.5;
        fig_height = 5;
    end
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    assert(~any(isnan(ceres)));
    
    c1 = [1 0 0;
          0 0 1;
          0 1 0;
          1 1 0;
          0 1 1;
          1 0 1];
    c2 = unique(lines, 'rows', 'stable');
    c3 = [128 128 128;
          233 185 210;
          240 185 140;
          163 181 141;
           83 177 247;
          169 145 145;
          212  94  94;
          152 174 246;
           91 165 100]/255;
    
    colors = [c2; c1; c3; 2*c1/3; 1-(1-c1)/3; [0 0 0]; [1 1 1]/3];

    sz = 1;
    
    [~, reorder] = sort(plot_order_bg_to_fg);    
    
    hold on;
    for i = plot_order_bg_to_fg'
        if(any(ceres==i))
            scatter(x(ceres==i), y(ceres==i), sz, colors(i,:), 'filled', 'DisplayName', ceres_names{i});        
        else
            scatter(NaN, NaN, sz, colors(i,:), 'filled', 'DisplayName', ceres_names{i});
        end
    end
    hold off;
    dat_series = get(gca, 'Children');
    dat_series = flipud(dat_series(1:length(reorder)));
    dat_series = dat_series(reorder);
    
    if (nargin > 8)
        dat_series = dat_series(legend_reorder);
    end
    
    axis tight;
    axis off;
    box off;
  
    top_margin = 0.3;
    bottom_margin = 0;    
    if (body)
        left_margin = 1.5;
        right_margin = 0;
    else
        left_margin = 0.1;
        right_margin = 3;
    end
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-right_margin-left_margin fig_height-top_margin-bottom_margin]);
    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal', 'interpreter', 'none');

    lgd = legend(dat_series, 'FontSize', 5);
    lgd.Units = 'centimeters';
    
    if (double_column)
        lgd.NumColumns = 2;    
        lgd_height = fig_height-bottom_margin;        
    else
        lgd.NumColumns = 1;
        lgd_height = fig_height/2-bottom_margin;
    end
    horz_offset = 0.2;
    if (body)
        lgd.Position = [-horz_offset, fig_height-lgd_height, left_margin, lgd_height];
    else        
        lgd.Position = [fig_width-right_margin-horz_offset, fig_height-lgd_height, right_margin, lgd_height];
    end
   
    legend('boxoff');

end