%% LT 10/8/15 - each point is paired vectors (X, Y). will plot in square and plot 45 degress line
function hplot=lt_plot_45degScatter(X, Y, color, PlotOn, normalplot)

if ~exist('normalplot', 'var')
    normalplot = 1
end
   
%% defualts

if ~exist('PlotOn','var');
    PlotOn=1;
end


if ~exist('color','var');
    color='k';
end


%% run
X=reshape(X,1, numel(X));
Y=reshape(Y,1, numel(Y));

lowerlim=min([X Y]);
upperlim=max([X Y]);

xlim([lowerlim upperlim]);
ylim([lowerlim upperlim]);

line([lowerlim upperlim], [lowerlim upperlim]);
line([0 0], ylim);
line(xlim, [0 0]);

% if plot
if PlotOn==1
    if normalplot==1
        plot(X,Y, 'o', 'Color', color);
    else
    lt_plot(X, Y, {'Color',color});
    end
end
