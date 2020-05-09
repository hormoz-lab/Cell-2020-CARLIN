function sp = plot_allele_membership_overlap(bank, Nr, Nc, which_sp)

    assert(isa(bank, 'Bank'));
    N_samples = size(bank.allele_breakdown_by_sample,2);
    assert(N_samples == 3);
  
    if (nargin == 1)
        fig_width  = 2.75;
        fig_height = 5;
        figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
               'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
        sp = subplot(1,1,1);
    else
        sp = subplot(Nr, Nc, which_sp);
    end
    
    R = 0.65;
    sep = 0.55;
    
    w = 2*R+sep;
    h = 2*R+sqrt(3)*sep/2;
    x_offset = (fig_width-w)/2;
    
    lgd_height = 0.9;
    lgd_width  = 0.8;
      
    label_width = 2.5;
    label_height = 0.1;
    tight_margin = 0.1;

    centre = [x_offset + w/2 - sep/2, lgd_height  + R + sqrt(3)/2*sep;              
              x_offset + w/2 + sep/2, lgd_height  + R + sqrt(3)/2*sep;
              x_offset + w/2        , lgd_height  + R;              
              x_offset + w/2 - sep/2, lgd_height + label_height + 2*tight_margin + h + R + sqrt(3)/2*sep;
              x_offset + w/2 + sep/2, lgd_height + label_height + 2*tight_margin + h + R + sqrt(3)/2*sep; 
              x_offset + w/2        , lgd_height + label_height + 2*tight_margin + h + R];
       
    theta = linspace(0, 2*pi, 360);
    x = R*cos(theta);
    y = R*sin(theta);
    
    hold on;
    patch(centre(1,1)+x, centre(1,2)+y, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(centre(2,1)+x, centre(2,2)+y, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(centre(3,1)+x, centre(3,2)+y, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(centre(4,1)+x, centre(4,2)+y, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(centre(5,1)+x, centre(5,2)+y, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(centre(6,1)+x, centre(6,2)+y, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold off;
    
    ax = gca;
    ax.Units = 'normalized';
    ax.Position = [0 0 1 1];
    ax.Color = 'none';
    ax.XColor = 'none';
    ax.YColor = 'none';    
    ax.Visible = 'off';
    xlim([0 fig_width]);
    ylim([0 fig_height]);
        
    sample_indicator = logical([1 0 0;
                                0 1 0;
                                0 0 1;
                                1 1 0;
                                1 0 1;
                                0 1 1;
                                1 1 1]);    
    [is, allele_membership] = ismember(logical(bank.allele_breakdown_by_sample), sample_indicator, 'rows');    
    assert(all(is));
    transcript_membership = arrayfun(@(i) sum(sum(bank.allele_breakdown_by_sample(allele_membership==i,:))), [1:7]');
    transcript_membership(7) = transcript_membership(7)-bank.summary.allele_freqs(1);
    allele_membership = accumarray(allele_membership,1);
    allele_membership(7) = allele_membership(7)-1;
    
    an_height = 0.1;
    an_width = 0.3;
    
    for i = 1:7
        an{i} = annotation('textbox', 'EdgeColor','none','string', num2str(transcript_membership(i)), ...
                           'Margin', 0, 'VerticalAlignment', 'middle', 'FitBoxToText', 'On', 'FontSize', 5);
        an{i}.Units = 'centimeters';
    end
    
    top_offset = 0.55;
    bottom_offset = 0.1;
    horz_offset = 0.1;
    
    an{1}.HorizontalAlignment = 'left';   an{1}.VerticalAlignment = 'top';    
    an{1}.Position = [x_offset+horz_offset, lgd_height+h-top_offset, an_width, an_height];
    
    an{2}.HorizontalAlignment = 'right';  an{2}.VerticalAlignment = 'top';    
    an{2}.Position = [x_offset+w-horz_offset, lgd_height+h-top_offset, an_width, an_height];
    
    an{3}.HorizontalAlignment = 'center'; an{3}.VerticalAlignment = 'bottom';
    an{3}.Position = [x_offset+w/2, lgd_height+bottom_offset, an_width, an_height];
    
    an{4}.HorizontalAlignment = 'center'; an{4}.VerticalAlignment = 'bottom';    
    an{4}.Position = (an{1}.Position+an{2}.Position)/2;    
    
    an{5}.HorizontalAlignment = 'right'; an{6}.VerticalAlignment = 'middle';
    an{5}.Position = (an{1}.Position+an{3}.Position)/2;
    
    an{6}.HorizontalAlignment = 'left'; an{6}.VerticalAlignment = 'middle';    
    an{6}.Position = (an{2}.Position+an{3}.Position)/2;
    
    an{7}.HorizontalAlignment = 'center'; an{7}.VerticalAlignment = 'top';    
    an{7}.Position = [sum(centre(1:3,:),1)/3-[0 an_height/2] an_width an_height];   
        
    an{5}.Position(1) = an{5}.Position(1) - 0.05;
    an{6}.Position(1) = an{6}.Position(1) + 0.05;
    an{4}.Position(2) = an{4}.Position(2) + 0.10;
    an{3}.Position(2) = an{3}.Position(2) + 0.10;
    
     for i = 1:7
        an{i+7} = annotation('textbox', 'EdgeColor','none','string', num2str(allele_membership(i)), ...
                           'Margin', 0, 'FitBoxToText', 'On', 'FontSize', 5);
        an{i+7}.Units = 'centimeters';
        an{i+7}.VerticalAlignment = an{i}.VerticalAlignment;
        an{i+7}.HorizontalAlignment = an{i}.HorizontalAlignment;
        an{i+7}.Position = an{i}.Position;
        an{i+7}.Position(2) = an{i}.Position(2) + h + label_height + 2*tight_margin;
     end
  
    allele_label = annotation('textbox', 'EdgeColor','none','string', 'Allele Membership', ...
                              'Margin', 0, 'FitBoxToText', 'On', 'FontSize', 6);
    allele_label.Units = 'centimeters';
    allele_label.HorizontalAlignment = 'center';
    allele_label.VerticalAlignment = 'bottom';    
    allele_label.Position = [x_offset+w/2-label_width/2 lgd_height+2*h+label_height+2*tight_margin label_width label_height];
    
    transcript_label = annotation('textbox', 'EdgeColor','none','string', 'Corresponding Transcripts', ...
                                   'Margin', 0, 'FitBoxToText', 'On', 'FontSize', 6);
    transcript_label.Units = 'centimeters';
    transcript_label.HorizontalAlignment = 'center';
    transcript_label.VerticalAlignment = 'middle';    
    transcript_label.Position = [x_offset+w/2-label_width/2 lgd_height+h+tight_margin label_width label_height];

    lgd = legend({'1'; '2'; '3'}, 'FontSize', 5, 'Location', 'North');
    title(lgd, 'Mouse', 'FontSize', 6, 'FontWeight', 'normal');
    lgd.Units = 'centimeters';
    lgd.NumColumns = 2;
    lgd.Position(1) = (fig_width-lgd_width)/2;
    lgd.Position(2) = 0;
    lgd.Position(3) = lgd_width;
    lgd.Position(4) = lgd_height;
    legend('boxoff');
    
    axis off;
    box off;
      
end