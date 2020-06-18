function [] = plot_rank(x, varargin)
loglog(sort(x, 'descend'), varargin{:});
xlabel('Rank');
ylabel('Frequency');
