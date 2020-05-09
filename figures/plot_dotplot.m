function plot_dotplot(dp_data, marker_genes, phenotype, title_str)

    louvain_cluster = horzcat(phenotype{:,2})';
    
    N_louv = cellfun(@length, phenotype(:,2));
    phenotype = phenotype(N_louv>0,1);
    N_louv = N_louv(N_louv>0);
    end_louv = cumsum(N_louv);
    start_louv = [1; end_louv(1:end-1)+1];
    end_louv = end_louv+0.25;
    start_louv = start_louv-0.25;
    
    fig_width  = 8.5;
    fig_height = 7;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);

    white_thresh = 0.9;
    cmap = [[linspace(0, white_thresh, 32) linspace(white_thresh, 1, 32)]
            [linspace(0, white_thresh, 32) linspace(white_thresh, 0, 32)];
            [linspace(1, white_thresh, 32) linspace(white_thresh, 0, 32)]]';
    colormap(cmap);
    sf = 10;
    
    hold on;
    
    plot(NaN, NaN, 'ko', 'markersize', sqrt((0+1)/sf), 'MarkerFaceColor', 'black');
    plot(NaN, NaN, 'ko', 'markersize', sqrt((25+1)/sf), 'MarkerFaceColor', 'black');
    plot(NaN, NaN, 'ko', 'markersize', sqrt((50+1)/sf), 'MarkerFaceColor', 'black');
    plot(NaN, NaN, 'ko', 'markersize', sqrt((75+1)/sf), 'MarkerFaceColor', 'black');
    plot(NaN, NaN, 'ko', 'markersize', sqrt((100+1)/sf), 'MarkerFaceColor', 'black');    
    scatter(dp_data.WhichMarkerGenes, dp_data.WhichLouvainCluster, (dp_data.PctExp+1)/10, dp_data.AvgExpScaled, 'filled');
    hold off;
    axis tight;
    xb = [1 length(marker_genes)];
    xb = [xb(1)-0.5 xb(2)+1.0];
    xlim([xb(1) xb(2)]);
    yb = ylim;
    yb = [yb(1)-0.5 yb(2)+1.0];
    ylim([yb(1) yb(2)]);  
    set(gca, 'Xtick', [1:length(marker_genes)]', 'Xticklabel', marker_genes);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xtickangle(90);
    xlabel('Marker Gene', 'FontSize', 6);
    set(gca, 'YTick', [1:length(louvain_cluster)], 'Yticklabel', cellstr(num2str(louvain_cluster)));
    set(get(gca, 'YAxis'), 'Dir', 'reverse', 'FontSize', 5);
    ylabel('Louvain Cluster', 'FontSize', 6);
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    caxis([-2.5 2.5]);
        
    yyaxis right;
    for i = 1:length(N_louv)
        line([xb(2) xb(2)], [start_louv(i) end_louv(i)], 'Color', 'black', 'LineWidth', 2.0);
    end    
    set(gca, 'Ytick', (start_louv+end_louv)/2, 'Yticklabel', phenotype, 'YColor', 'black');
    yrulers = get(gca, 'YAxis');        
    set(yrulers(2), 'Dir', 'reverse', 'FontSize', 5);
    set(gca, 'TickDir', 'out', 'LineWidth', 0.1);
    ylabel('Cell Type', 'FontSize', 6);
    ylim([yb(1) yb(2)]);
    
    cb = colorbar('Location', 'EastOutside', 'FontSize', 5);
    cb.AxisLocation = 'out';
    title(cb, 'AvgExpScaled', 'FontSize', 6, 'FontWeight', 'normal');
    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    box off;
    
    left_margin = 0.8;
    top_margin = 0.5;
    bottom_margin = 1.0;
    right_margin = 1.0;
    legend_width = 1.2;
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, ...
                                                 fig_width-left_margin-right_margin-legend_width fig_height-top_margin-bottom_margin]);

    leg_height = 1.5;
    scale_width = 0.3;
    tight_margin = 0.1;
    cb_height = fig_height-top_margin-leg_height-2*tight_margin;
                                             
    cb.Units = 'centimeters';
    cb_pos = cb.Position;
    cb.Position(2) = leg_height+2*tight_margin;
    cb.Position(3) = scale_width;
    cb.Position(4) = cb_height;
    cb.TickLabels(end) = {['>= ' cb.TickLabels{end}]};
    cb.TickLabels(1)   = {['<= ' cb.TickLabels{1}]};
    
    lgd = legend({'0%'; '25%'; '50%'; '75%'; '100%'}, 'Location', 'SouthEastOutside', 'FontSize', 5, 'Color', 'black');
    lgd.Units = 'centimeters';    
    lgd.Position(1) = cb.Position(1)-tight_margin;
    lgd.Position(2) = tight_margin;
    lgd.Position(3) = cb.Position(3);
    lgd.Position(4) = leg_height;
    
    title(lgd, '% Expressed', 'FontSize', 6', 'FontWeight', 'normal');
    legend boxoff;
    
end
