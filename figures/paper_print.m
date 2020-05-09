function paper_print(filename, opts)

    if (nargin == 1)     
        try
            print('-dtiff','-r4800', sprintf('%s.tiff', filename));
        catch
            print('-dtiff','-r2400', sprintf('%s.tiff', filename));
        end        
    elseif (strcmp(opts, 'SVG'))
        print(sprintf('%s.svg', filename), '-dsvg');
    elseif (strcmp(opts, 'PNG'))
        print(sprintf('%s.png', filename), '-dpng', '-r2400');
    end
    
end
    
        