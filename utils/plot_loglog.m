function [h2, centerbins, freq] = plot_loglog(x, varargin)

linespec = 'o';
step = 1;
xmin = min(x);
xmax = max(x);
show = true;
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'linespec'
            linespec = varargin{2};
        case 'step'
            step = varargin{2};
        case 'xmin'
            xmin = varargin{2};
        case 'xmax'
            xmax = varargin{2};
        case 'show'
            show = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}]);
    end
    varargin(1:2) = [];
end

edgebins = 2.^(log2(xmin):step:log2(xmax));
sizebins = edgebins(2:end) - edgebins(1:end-1);
sizebins(end+1) = 1;
centerbins = edgebins;
counts = histc(x, edgebins);
freq = counts./sizebins/length(x);
if show
    h2 = loglog(centerbins, freq, linespec, 'linewidth', 1.5);
    xlabel('Size', 'fontsize', 16)
    ylabel('Distribution', 'fontsize', 16)
else
    h2 = 0;
end
