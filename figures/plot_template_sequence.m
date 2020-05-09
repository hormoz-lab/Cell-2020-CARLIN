function plot_template_sequence(which)

    ref = CARLIN_def.getInstance;        
    [~, template] = CARLIN_def.cas9_align(ref.seq.CARLIN);
    RGB = get_sequence_coloring({template}, 'bp');
    
    fig_height = 0.2;
    fig_width = 10.8;    
    
    if (strcmp(which, 'All'))
        N_rep = [10 1];
        bounds = [ref.bounds.prefix(1) ref.bounds.postfix(2)];
    else
        N_rep = [10 5];
        if (strcmp(which, 'Repeat'))            
            bounds = [ref.bounds.consites(1,1) ref.getInstance.bounds.pams(1,2)];
        elseif (strcmp(which, 'Prefix'))            
            bounds = [ref.bounds.prefix(1) ref.getInstance.bounds.prefix(2)];
        elseif (strcmp(which, 'Postfix'))            
            bounds = [ref.bounds.postfix(1) ref.getInstance.bounds.postfix(2)];
        end
    end
    
    fig_width = fig_width*(2+(diff(bounds)+1)*N_rep(2))/(2+ref.width.CARLIN);
        
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);    
    sp = subplot(1, 1, 1);
    set(sp, 'Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height]);       
    
    imshow(padarray(repelem(RGB(:,bounds(1):bounds(2),:), N_rep(1), N_rep(2), 1), [1 1], 0)); hold on;
    h = imshow(padarray(repmat(repelem(ref.alpha.CARLIN(bounds(1):bounds(2)), N_rep(2)), [N_rep(1), 1]), [1 1], 0));  hold off;
    if (~strcmp(which, 'All'))
        set(h, 'AlphaData', 0.8);
    end

    box on;
    daspect([(fig_height)/(N_rep(1)+2), (fig_width)/(2+(diff(bounds)+1)*N_rep(2)), 1]);

end