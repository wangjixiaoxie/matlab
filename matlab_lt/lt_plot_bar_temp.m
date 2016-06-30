%% 7/21/15 - LT plots bar() but with transparent face, color, and width 0.5
function hbar=lt_plot_bar(X, Y, varargin)

% E.g.
% lt_bar(X, Y, {'Color', 'k'}):
% Other modifiers;
% {'BarWidth', 0.5}; (default 0.5)
% {'Errors', Yerr}, is for errorbars. Yerr is vector with errors
%   (symmetrical) NOTE: if using errorbars, must supply an X, but X will
%   not be used (will plot Y in order)

%%
    % default
    color='w';
    barwidth=0.6;
    ploterror=0;

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
    
    % error
    ind=find(strcmp(varargin{1},'Errors'));
    if isempty(ind);
        ploterror=0;
    else
        ploterror=1;
        Yerr=varargin{1}{ind+1};
    end
end

if ploterror==0;
% 1) no errorbar
hbar=bar(X, Y, barwidth);
set(hbar, 'FaceColor',color)
% set(hbar, 'LineWidth', 2)
else
% 2) with errorbar
hbar=barwitherr(Yerr, Y);
set(hbar, 'FaceColor',color);
set(hbar, 'BarWidth', barwidth);

end


% make translucent
hbar_ch=get(hbar, 'Children');
set(hbar_ch, 'FaceAlpha', 0.3);
