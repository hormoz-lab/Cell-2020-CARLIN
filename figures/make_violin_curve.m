function [xfine, yfine] = make_violin_curve(x, y, which)

    if (nargin == 2)
        which = 'P';
    end
    assert(issorted(y));

    y = [y; y(2:end)-0.5; y+0.5];
    x = padarray(x, length(y)-length(x), 0, 'post');
    [~, idx] = sort(y, 'ascend');
    y = y(idx); x = x(idx);
    u = unique([y x], 'rows', 'stable');
    y = u(:,1); x = u(:,2);
    
    yfine = [];
    xfine = [];
    N_xfine = 10;
    for i = 1:length(x)-1
        if (x(i) ~= x(i+1))
            xmid = (x(i)+x(i+1))/2;
            ymid = (y(i)+y(i+1))/2;
            if (which == 'P')
                if (x(i) > 0 && x(i+1) == 0)
                    yvec = [x(i+1); xmid; x(i); 0];
                    xmat = [y(i+1)^3   y(i+1)^2 y(i+1) 1; 
                            ymid^3     ymid^2   ymid   1;
                            y(i)^3     y(i)^2   y(i)   1;
                            3*y(i+1)^2 2*y(i+1) 1      0];
                    coeffs = inv(xmat)*yvec;
                    yy = linspace(y(i+1), y(i), N_xfine)';
                    xx = [yy.^3 yy.^2 yy.^1 yy.^0]*coeffs;                
                else
                    yvec = [x(i); xmid; x(i+1); 0];
                    xmat = [y(i)^3   y(i)^2   y(i)   1; 
                            ymid^3   ymid^2   ymid   1;
                            y(i+1)^3 y(i+1)^2 y(i+1) 1;
                            3*y(i)^2 2*y(i)   1      0];
                    coeffs = inv(xmat)*yvec;
                    yy = linspace(y(i), y(i+1), N_xfine)';
                    xx = [yy.^3 yy.^2 yy.^1 yy.^0]*coeffs;
                end
                xfine = [xfine; xx];
                yfine = [yfine; yy];
            elseif (which == 'C')
                a = abs(xmid-x(i));
                b = abs(ymid-y(i));
                xx = linspace(0, a, N_xfine)';
                yy = sqrt((1-xx.^2/a^2)*b^2);
                if (x(i) > 0 && x(i+1) == 0)
                    xfine = [xfine; xmid+xx; xmid-xx];
                    yfine = [yfine; y(i)+yy; y(i+1)-yy];
                elseif (x(i) == 0 && x(i+1) > 0)
                    xfine = [xfine; xmid-xx; xmid+xx];
                    yfine = [yfine; y(i)+yy; y(i+1)-yy];
                end
            end         
        end
    end
    u = unique([yfine xfine], 'rows');
    yfine = u(:,1); 
    xfine = u(:,2);
    [~, idx] = sort(yfine, 'ascend');
    yfine = yfine(idx);
    xfine = xfine(idx);
    xfine = xfine/2;
    xfine = [xfine; -flipud(xfine)];
    yfine = [yfine;  flipud(yfine)];
    
end