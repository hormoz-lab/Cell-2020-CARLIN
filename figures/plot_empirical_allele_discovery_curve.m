function sp = plot_empirical_allele_discovery_curve(bank, Nr, Nc, which_sp)

    assert(isa(bank, 'Bank'));

    p = repelem([1:size(bank.summary.alleles,1)-1]', bank.summary.allele_freqs(2:end));
    
    Npoints_fine = 101;
    Ntrials = 100;
    
    n_fine = floor(linspace(100, bank.model.transcripts.edited, Npoints_fine));
    n_diversity = zeros(Ntrials, Npoints_fine);
    
    for i = 1:Ntrials
        for j = 1:Npoints_fine
            datavector = datasample(p, n_fine(j), 'Replace', false);
            n_diversity(i,j) = length(unique(datavector));
        end
    end
    
    Npoints_coarse = 11;
    Ninterp = 10;
    
    n_coarse = floor(linspace(100, bank.model.transcripts.edited, Npoints_coarse));
    m = arrayfun(@(n) linspace(1, n*log(n)-1, Ninterp)', n_coarse, 'un', false);
    m = horzcat(m{:});
    m = floor(m);
    
    n_extrapolated = zeros(Ninterp, Npoints_coarse, Ntrials);
    
    for i = 1:Ninterp
        for j = 1:Npoints_coarse
            for k = 1:Ntrials
                datavector = datasample(p, n_coarse(j), 'Replace', false);
                freq_counts = accumarray(nonzeros(accumarray(datavector,1)),1);                
                n_extrapolated(i,j,k) = Orlitzky(m(i,j), freq_counts);
            end
        end
    end
        
    if (nargin == 1)
        fig_width  = 4.275;
        fig_height = 5;
        figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
               'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
        sp = subplot(1,1,1);
    else
        sp = subplot(Nr, Nc, which_sp);
    end
    
    cmap = colormap('jet');
    cmap = cmap(1:5:end,:);
    cmap = cmap(3:end,:);
    
    hold on;
    patch([0 0 bank.model.transcripts.edited bank.model.transcripts.edited], [0 bank.model.alleles.observed bank.model.alleles.observed 0], ...
          'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none');    
    plot(n_fine, mean(n_diversity,1), 'LineWidth', 1, 'Color', 'black');
    matching = find(ismember(n_fine, n_coarse));
    for i = 1:Npoints_coarse        
        errorbar(m(:,i)+n_coarse(i), mean(n_extrapolated(:,i,:),3)+mean(n_diversity(:,matching(i))), ...
            std(n_extrapolated(:,i,:),[],3), 'color', cmap(i,:), 'LineStyle', '--');
        scatter(n_fine(matching(i)), mean(n_diversity(:,matching(i))), 20, cmap(i,:), 'filled');
    end
    hold off;
    
    axis tight;
    
    max_expo = floor(log10(bank.model.alleles.CI_95_UB));  
    set(get(gca, 'YAxis'), 'Exponent', max_expo, 'FontSize', 5);    
    ylabel('Alleles', 'FontSize', 6);
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('Edited Transcripts', 'FontSize', 6);    
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    box off;
    
    title('Non-Parametric', 'FontSize', 6, 'FontWeight', 'normal');    
    
    left_margin = 0.6;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    tight_margin = 0.1;        
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-tight_margin fig_height-top_margin-bottom_margin]);
    
    lgd = legend({'Observed'; 'Interpolated'; sprintf('Extrapolated\n(Mean +/- SD')}, 'Location', 'Southeast', 'FontSize', 5);
    lgd.Units = 'centimeters';        
    leg_width = 3.0;
    leg_height = 1.2;
    lgd.Position = [fig_width-leg_width, bottom_margin+5*tight_margin, leg_width, leg_height];
    
    legend('boxoff');
    
end
