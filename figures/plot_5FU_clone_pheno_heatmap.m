function plot_5FU_clone_pheno_heatmap(samples, coarse_grain_pheno, title_str, colors, N_downsample)

    % Expects samples{1} to be control which needs downsampling and
    % samples{2} to be +5FU

    fig_width = 6.0;
    fig_height = 12;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
    
    N_mice = 2;
    N_tissues = length(coarse_grain_pheno);
    assert(length(title_str) == N_mice);
    assert(length(samples)   == N_mice);
       
    top_margin = 0.3;
    bottom_margin = 0.7;
    right_margin = 0.4;
    left_margin = 1.1;
   
    cb_width = 0.2;
    tight_margin = 0.2;
    
    heatmap_width = (fig_width-N_mice*left_margin-right_margin-cb_width-tight_margin)/2;
    rel_height = fig_height-top_margin-bottom_margin;
    
    sp = cell(N_mice,1);
    
    cmap = ones(64,3);
    cmap(:,2) = linspace(0, 1, 64);
    cmap(:,3) = linspace(0, 1, 64);
    colormap(flipud(cmap));
    max_count = max(cellfun(@(x) max(max(x.clone_by_cg_pheno_breakdown(x.sig_alleles.all,:))), samples));
    
    if (nargin == 4)
        N_downsample = length(samples{2}.sig_alleles.all);
    end
    
    for i = 1:N_mice
        
        sp{i} = subplot(1,N_mice,i);
        
        [~, idx1] = sortrows(samples{i}.clone_by_cg_pheno_breakdown(samples{i}.sig_alleles.hsc_derived_allele,:), 'descend');
        [~, idx2] = sortrows(samples{i}.clone_by_cg_pheno_breakdown(samples{i}.sig_alleles.hsc_only_allele,:), 'descend');
        [~, idx3] = sortrows(samples{i}.clone_by_cg_pheno_breakdown(samples{i}.sig_alleles.derived_only_allele,:), 'descend');
        
        N_alleles = length(samples{i}.sig_alleles.all);
        reordered_alleles = [samples{i}.sig_alleles.hsc_derived_allele(idx1);
                             samples{i}.sig_alleles.hsc_only_allele(idx2); 
                             samples{i}.sig_alleles.derived_only_allele(idx3)];
        assert(length(reordered_alleles) == N_alleles);
        
        if (i == 1 || nargin == 5)
            % Subsample in a way that preserves matrix structure
        
            rng(78436743);
            ss_factor = N_downsample / length(samples{i}.sig_alleles.all);
            
            pheno_indicator = logical(dec2bin([1:2^length(coarse_grain_pheno)-1])-'0');
            [is, where] = ismember(logical(samples{i}.clone_by_cg_pheno_breakdown(reordered_alleles,:)), pheno_indicator, 'row');
            ss = arrayfun(@(i) find(where==i), [1:size(pheno_indicator,1)]', 'un', false);
            ss = cellfun(@(x) datasample(x, round(length(x)*ss_factor), 'Replace', false), ss, 'un', false);
            ss = vertcat(ss{:});
            % If we can't match the exact length, take from the top (HSC)
            % or remove from the bottom (Other), since these are most
            % interesting.
            if (length(ss) > N_downsample)
                ss = ss(1:N_downsample);
            else
                pad = setdiff([1:length(reordered_alleles)]', ss);
                ss = sort([ss; pad(1:(N_downsample-length(ss)))]);
            end
            reordered_alleles = reordered_alleles(ss);
            N_alleles = length(reordered_alleles);
            assert(N_alleles == N_downsample);
        end 
        
        imagesc(samples{i}.clone_by_cg_pheno_breakdown(reordered_alleles,:), 'XData', 1:N_tissues, 'YData', 1:N_alleles);        
        caxis([0 max_count]);
        set(gca, 'Xtick', 1:N_tissues, 'xticklabel', coarse_grain_pheno, 'TickDir', 'out');
        set(gca, 'YTick', 1:N_alleles, 'yticklabel', arrayfun(@(x) sprintf('Clone %-3d', x), reordered_alleles, 'un', false));
        ax = gca;
        
        for j = 1:length(ax.YTickLabel)            
            if (ismember(reordered_alleles(j), samples{i}.sig_alleles.hsc_derived_allele))
                color_prefix = sprintf('\\color[rgb]{%f,%f,%f}', colors(1,1), colors(1,2), colors(1,3));                
            elseif (ismember(reordered_alleles(j), samples{i}.sig_alleles.hsc_only_allele))                
                color_prefix = sprintf('\\color[rgb]{%f,%f,%f}', colors(2,1), colors(2,2), colors(2,3));
            elseif (ismember(reordered_alleles(j), samples{i}.sig_alleles.derived_only_allele))
                color_prefix = sprintf('\\color[rgb]{%f,%f,%f}', colors(3,1), colors(3,2), colors(3,3));
            end        
            ax.YTickLabel{j} = sprintf('%s%s', color_prefix, ax.YTickLabel{j});
        end
        
        set(get(ax, 'XAxis'), 'FontSize', 5);
        set(get(ax, 'YAxis'), 'FontSize', 5);
        xtickangle(90);
        title(title_str{i}, 'FontSize', 6, 'FontWeight', 'normal');
        bg{i} = gca;    
        set(bg{i},'box','off','color','none')        
        fg{i} = axes('Position', get(bg{i},'Position'),'box','on','xtick',[],'ytick',[], 'LineWidth', 1.0);    
        axes(bg{i});
    end
    
    cb = colorbar('Location', 'EastOutside', 'FontSize', 5);
    cb.AxisLocation = 'out';
    
    set(sp{1}, 'Units', 'centimeters', 'Position', [  left_margin bottom_margin heatmap_width rel_height]);
    set(fg{1}, 'Units', 'centimeters', 'Position', [  left_margin bottom_margin heatmap_width rel_height]);
    set(sp{2}, 'Units', 'centimeters', 'Position', [2*left_margin+heatmap_width bottom_margin heatmap_width rel_height]);
    set(fg{2}, 'Units', 'centimeters', 'Position', [2*left_margin+heatmap_width bottom_margin heatmap_width rel_height]);
    
    cb.Units = 'centimeters';
    cb.Position(1) = 2*left_margin+2*heatmap_width+tight_margin;
    cb.Position(2) = bottom_margin;
    cb.Position(3) = cb_width;
    cb.Position(4) = rel_height;
    title(cb, 'Cells', 'FontSize', 6, 'FontWeight', 'normal', 'VerticalAlignment', 'baseline');
    
end   