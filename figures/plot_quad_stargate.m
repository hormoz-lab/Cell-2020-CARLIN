function fig = plot_quad_stargate(summary, labels)

    assert(isa(summary{1}, 'ExperimentSummary') && size(summary,1)==4);
    
    chord_radius = 3;
    inner_margin = 0.2;
    chord_label_height = 0.2;
    chord_label_width = 0.3;
    
    scale_width = 0.7;
    scale_label_width = 1.0;
    scale_label_height = 0.8;
        
    fig_width  = 2*chord_radius + 2*inner_margin + scale_width;
    fig_height = 2*chord_radius + 1*inner_margin;
    
    fig = figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
                 'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
    
    chord_sp = arrayfun(@(i) plot_stargate.create(summary{i}, 2, 2, i), [1:4]', 'un', false);
    
    linkaxes([chord_sp{1} chord_sp{3}], 'x');
    linkaxes([chord_sp{2} chord_sp{4}], 'x');
    linkaxes([chord_sp{1} chord_sp{2}], 'y');
    linkaxes([chord_sp{3} chord_sp{4}], 'y');
    
    set(chord_sp{1}, 'Units', 'centimeters', 'Position', [0*chord_radius + 0*inner_margin 1*chord_radius + 1*inner_margin chord_radius chord_radius]);
    set(chord_sp{2}, 'Units', 'centimeters', 'Position', [1*chord_radius + 1*inner_margin 1*chord_radius + 1*inner_margin chord_radius chord_radius]);
    set(chord_sp{3}, 'Units', 'centimeters', 'Position', [0*chord_radius + 0*inner_margin 0*chord_radius + 0*inner_margin chord_radius chord_radius]);
    set(chord_sp{4}, 'Units', 'centimeters', 'Position', [1*chord_radius + 1*inner_margin 0*chord_radius + 0*inner_margin chord_radius chord_radius]);

    if (nargin == 2)
        assert(size(labels,1)==4)
        sg_an = cell(4,1);
        for i = 1:4
            sg_an{i} = annotation('textbox', 'EdgeColor', 'none', 'string', labels{i}, 'FitBoxToText', 'On', 'Margin', 1, 'FontSize', 6);
            sg_an{i}.Units = 'centimeters';
        end
        sg_an{1}.Position = [0*chord_radius + 0*inner_margin fig_height-chord_label_height   chord_label_width chord_label_height];
        sg_an{2}.Position = [1*chord_radius + 1*inner_margin fig_height-chord_label_height   chord_label_width chord_label_height];
        sg_an{3}.Position = [0*chord_radius + 0*inner_margin chord_radius-chord_label_height chord_label_width chord_label_height];
        sg_an{4}.Position = [1*chord_radius + 1*inner_margin chord_radius-chord_label_height chord_label_width chord_label_height];        
    end

    range_line_width_pts = (plot_stargate.max_width_in_pts - 1 - plot_stargate.min_width_in_pts);
    orders_of_magnitude = floor(range_line_width_pts)+1;
    
    points_to_cm = double(unitConversionFactor("point", "cm"));
    
    % Max line width defined in plot_stargate corresponds to 0 = log_2(x=1)
    % We set an upper bound for the legend of x=0.5 => log_2(0.5)=-1
    max_line_width = (plot_stargate.max_width_in_pts-1) * points_to_cm;
    min_line_width = plot_stargate.min_width_in_pts * points_to_cm;
    
    assert(max_line_width <= scale_width);
    
    scale_height = fig_height - inner_margin - 2*scale_label_height;
    
    an = annotation('textbox', 'Edgecolor', 'none', 'string', 'log_2 Cell Fraction', ...
                    'HorizontalAlignment', 'right', 'FitBoxToText', 'On', 'Margin', 1, 'Fontsize', 6);
    an.Units = 'centimeters';
    an.Position = [fig_width-inner_margin-scale_label_width, fig_height-scale_label_height, scale_label_width, scale_label_height];
    
    ax = axes('Units', 'centimeters', 'Position', [2*chord_radius + 2*inner_margin, inner_margin+scale_label_height, max_line_width, scale_height]);
    patch([max_line_width-min_line_width, 0, max_line_width, max_line_width], ...
          [0 scale_height scale_height 0], [0.5 0 0.5], 'EdgeColor', 'none');
    set(gca, 'YAxisLocation', 'right');
    set(get(gca, 'YAxis'), 'FontSize', 5);
    
    yticks(scale_height-(scale_height/range_line_width_pts)*[orders_of_magnitude-1:-1:0]);
    yticklabels(cellstr(num2str([-orders_of_magnitude:1:-1]')));
    axis tight; box off; 
    ax.XAxis.Visible = 'off';
    
    legend_label_height = 0.15;
    legend_label_width = 0.9;
    legend_box_width = 0.2;
    legend_horz_spacer = 0.05;
    legend_vert_spacer = 0.05;
    
    an = annotation('rectangle', 'Edgecolor', 'b', 'FaceColor', 'b');
    an.Units = 'centimeters';
    an.Position = [fig_width-legend_horz_spacer-legend_box_width-legend_label_width, legend_label_height+inner_margin+legend_vert_spacer, legend_box_width, legend_label_height];
    
    an = annotation('textbox', 'Edgecolor', 'none', 'string', 'Insertion', ...
                    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'bottom', 'FitBoxToText', 'On', 'Margin', 0, 'Fontsize', 6);
    an.Units = 'centimeters';
    an.Position = [fig_width-legend_label_width, legend_label_height+inner_margin, legend_label_width, legend_label_height];
    
    
    an = annotation('rectangle', 'Edgecolor', 'r', 'FaceColor', 'r');
    an.Units = 'centimeters';
    an.Position = [fig_width-legend_horz_spacer-legend_box_width-legend_label_width, legend_vert_spacer, legend_box_width, legend_label_height];
    
    
    an = annotation('textbox', 'Edgecolor', 'none', 'string', 'Deletion', ...
                    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'bottom', 'FitBoxToText', 'On', 'Margin', 0, 'Fontsize', 6);
    an.Units = 'centimeters';
    an.Position = [fig_width-legend_label_width, 0, legend_label_width, legend_label_height];
    
end
