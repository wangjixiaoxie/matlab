%% LT 10/8/15 - each point is paired vectors (X, Y). will plot in square and plot 45 degress line
function lt_plot_makesquare_plot45line(haxes, color)

%% defualts


if ~exist('color','var');
    color='k';
end


%% run

Ylim=get(haxes, 'Ylim');
Xlim=get(haxes, 'Xlim');


lowerlim=min([Ylim(1) Xlim(1)]);
upperlim=max([Ylim(2) Xlim(2)]);

amount_to_add=(upperlim-lowerlim)/10;
lowerlim=lowerlim+amount_to_add;
upperlim=upperlim+amount_to_add;


xlim([lowerlim upperlim]);
ylim([lowerlim upperlim]);

line([lowerlim upperlim], [lowerlim upperlim], 'Color',color);
line([0 0], ylim);
line(xlim, [0 0]);
