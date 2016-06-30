%% 7/21/15 - LT plots bar() but with transparent face, color, and width 0.5
function hbar=lt_plot_bar(X, Y, varargin)

% E.g.
% lt_bar(X, Y, {'Color', 'k'}):
% Other modifiers;
% {'BarWidth', 0.5}; (default 0.5)
% {'Errors', Yerr}, is for errorbars. Yerr is vector with errors
%   (symmetrical) NOTE: if using errorbars, must supply an X, if X is
%   length 1 then will ignore X and just plot at X=1 (one time I tried using X when X was length 1 and
%   was weird)
% {'FaceAlpha', 0.3}  - transparency

%%
    % default
    color='k';
    barwidth=0.6;
    ploterror=0;
    FaceAlpha=0.2;
    
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
    
    % transparency
    ind=find(strcmp(varargin{1},'FaceAlpha'));
    if ~isempty(ind);
        FaceAlpha=varargin{1}{ind+1};
    end
    
end

if ploterror==0;
    % 1) no errorbar
    hbar=bar(X, Y, barwidth);

else
    % 2) with errorbar
    if length(X)==1;
%         X=[X X+1];
%         Y=[Y nan];
%         Yerr=[Yerr 0];
        hbar=barwitherr_lt(Yerr, X, Y);
    else
        hbar=barwitherr_lt(Yerr, X, Y);
    end
    
    set(hbar, 'BarWidth', barwidth);
    
end

    set(hbar, 'FaceColor',color)
set(hbar, 'LineWidth', 1.5);

% make translucent
hbar_ch=get(hbar, 'Children');
set(hbar_ch, 'FaceAlpha', FaceAlpha);
hold on;
