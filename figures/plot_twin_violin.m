function plot_twin_violin(x1, y1, x2, y2, xlabels, ann_text)

    [xv1, yv1] = make_violin_curve(x1, y1, 'P');
    [xv2, yv2] = make_violin_curve(x2, y2, 'P');
    
    fig_width  = 3.7;
    fig_height = 3.7;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
       
    x_offset = max((max(xv1)+max(xv2))*1.2, 10);
    
    hold on;
    scatter(zeros(size(y1)), y1, 1, 'b', 'filled');
    scatter(x_offset*ones(size(y2)), y2, 1, 'b', 'filled');
    fill(xv1, yv1, 'b', 'EdgeColor', 'b', 'LineWidth', 0.5, 'FaceAlpha', 0.5);
    fill(xv2+x_offset, yv2, 'b', 'EdgeColor', 'b', 'LineWidth', 0.5, 'FaceAlpha', 0.5);
              
    axis tight;
    
    an_height = diff(ylim)*0.05;
    an_width = 600;
    
    y_offset = max(diff(ylim)*0.05, 2);
    maxy1 = max(y1)+y_offset;
    maxy2 = max(y2)+y_offset;
    maxy = max(maxy1, maxy2)+y_offset;
    
    line([0 0], [maxy1 maxy], 'Color', 'black', 'LineWidth', 0.5);
    line([x_offset x_offset], [maxy2 maxy], 'Color', 'black', 'LineWidth', 0.5);
    line([0 x_offset], [maxy maxy], 'Color', 'black', 'LineWidth', 0.5);
    line([1 1]*(x_offset)/2, [maxy maxy+y_offset], 'Color', 'black', 'LineWidth', 0.5);
    
    hold off;
    
    ylim([0 max(ylim)+an_height+y_offset]);
    xlim([min(xv1)-5 max(xv2)+x_offset+5]);
    
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    ylabel('Number of Cells', 'FontSize', 6);
        
    set(get(gca, 'XAxis'), 'FontSize', 6);
    set(gca, 'TickLabelInterpreter', 'tex')
    set(gca, 'Xtick', [0 x_offset], 'xticklabel', xlabels);
    
    title('HSC-Rooted Clone Sizes', 'FontWeight', 'normal', 'FontSize', 6);
    
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
   
    box off;
    
    left_margin = 0.7;
    top_margin = 0.4;
    bottom_margin = 0.6;
    right_margin = 0.3;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin]);
    
    pos = get(gca, 'Position');
    
    an = annotation('textbox', 'EdgeColor','none','string', ann_text, 'Margin', 0, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 5);    
    an.Position = [((x_offset/2-an_width/2 + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1))/fig_width,...      
                   ((maxy+y_offset - min(ylim))/diff(ylim) * pos(4) + pos(2))/fig_height, ...
                   ((an_width)/diff(xlim) * pos(3))/fig_width, ...
                   ((an_height)/diff(ylim) * pos(4))/fig_height];
end

