function sp = plot_theoretical_allele_discovery_curve(bank, Nr, Nc, which_sp)

    assert(isa(bank, 'Bank'));

    interp.transcripts = linspace(0, bank.model.transcripts.edited, 101);
    [interp.alleles.mu, interp.alleles.sig] = bank.interpolate_alleles(interp.transcripts);

    extrap.max_transcripts = bank.extrapolate_transcripts(0.99*bank.model.alleles.estimated);
    extrap.transcripts = linspace(bank.model.transcripts.edited, extrap.max_transcripts, 1000);
    [extrap.alleles.mu, extrap.alleles.CI_95_LB, extrap.alleles.CI_95_UB] = bank.extrapolate_alleles(extrap.transcripts);

    if (nargin == 1)
        fig_width  = 4.275;
        fig_height = 5;
        figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
               'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
        sp = subplot(1,1,1);
    else
        sp = subplot(Nr, Nc, which_sp);
    end
    
    hold on;
    patch([0 0 bank.model.transcripts.edited bank.model.transcripts.edited], [0 bank.model.alleles.observed bank.model.alleles.observed 0], ...
          'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    plot(interp.transcripts, interp.alleles.mu, 'color', 'black', 'LineWidth', 1.0); hold on;
    plot(extrap.transcripts, extrap.alleles.mu, 'color', 'red', 'LineWidth', 1.0);
    plot([0 extrap.transcripts(end)], [bank.model.alleles.estimated bank.model.alleles.estimated], ...
        'color', 'blue', 'LineWidth', 1.0, 'LineStyle', '--');
    patch([0 0 extrap.transcripts(end) extrap.transcripts(end)], ...
          [bank.model.alleles.CI_95_LB bank.model.alleles.CI_95_UB bank.model.alleles.CI_95_UB bank.model.alleles.CI_95_LB],  ...
          'blue', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    fill([extrap.transcripts fliplr(extrap.transcripts)], ...
         [extrap.alleles.CI_95_LB fliplr(extrap.alleles.CI_95_UB)], 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    
    hold off;
    axis tight;
    
    max_expo = floor(log10(bank.model.alleles.CI_95_UB));
    set(get(gca, 'YAxis'), 'Exponent', max_expo, 'FontSize', 5);
    ylim([0 ceil(bank.model.alleles.CI_95_UB/10^max_expo)*10^max_expo]);
    ylabel('Alleles', 'FontSize', 6);
    
    max_expo = floor(log10(extrap.transcripts(end)));
    set(get(gca, 'XAxis'), 'Exponent', max_expo, 'FontSize', 5);
    xlabel('Edited Transcripts', 'FontSize', 6);
        
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    title('Parametric', 'FontSize', 6, 'FontWeight', 'normal');
    
    box off;
    
    left_margin = 0.6;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    tight_margin = 0.1;        
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-tight_margin fig_height-top_margin-bottom_margin]);
    
    lgd = legend({'Observed'; 'Interpolated'; sprintf('Extrapolated\n(95%% CI)'); sprintf('Model Estimate\n(95%% CI)')}, ...
                 'Location', 'Southeast', 'FontSize', 5);
    lgd.Units = 'centimeters';        
    leg_width = 3.0;
    leg_height = 1.5;
    lgd.Position = [fig_width-leg_width, bottom_margin+5*tight_margin, leg_width, leg_height];
    
    legend('boxoff');
    
end