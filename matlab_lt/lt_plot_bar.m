%% 7/21/15 - LT plots bar() but with transparent face, color, and width 0.5
function hbar=lt_plot_bar(X, Y, varargin)

% E.g.
% lt_bar(X, Y, {'Color', 'k'}):
% Other modifiers;
% {'BarWidth', 0.5}; (default 0.5

%%
    % default
    color='k';
    barwidth=0.6;

if nargin>2;
    % color
    ind=find(strcmp(varargin{1}, 'Color'));
    if ~isempty(ind)
    color=varargin{1}{ind+1};
    end
    
    % width
    ind=find(strcmp(varargin{1}, 'BarWidth'));
    if ~isempty(ind)
    barwidth=varargin{1}{ind+1};
    end
end


hbar=bar(X, Y, barwidth);
set(hbar, 'FaceColor',color)
% set(hbar, 'LineWidth', 2)


% make translucent
hbar_ch=get(hbar, 'Children');
set(hbar_ch, 'FaceAlpha', 0.3);
