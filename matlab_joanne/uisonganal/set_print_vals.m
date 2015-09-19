function set_print_vals(fig_handle, orient, papertype, x_per_paperinch, y_per_paperinch)

%sets values of the selected figure and its current axis for printing
% orient = paper orientation
% papertype = paper type
% x_per_paperinch = number of units along the x axis per inch of printed figure
% y_per_paperinch = number of units along the y axis per inch of printed figure

%for now, set default margin sizes
left_margin = 1.0;
top_margin = 1.0;

%get axis that is current within the selected figure
axis_handle = get(fig_handle,'CurrentAxes');

%save unit values and set to inches
t_paperunits = get(fig_handle,'PaperUnits');
set(fig_handle,'PaperUnits','inches');
t_fig_units = get(fig_handle,'Units');
set(fig_handle,'Units','inches');
t_axis_units = get(axis_handle,'Units');
set(axis_handle,'Units','inches');

%get other printing parameters
paperposition = get(fig_handle,'PaperPosition');
fig_position = get(fig_handle,'Position'); 
fig_width = fig_position(3);
fig_paperwidth = paperposition(3);
fig_height = fig_position(4);
fig_paperheight = paperposition(4);

xlim = get(axis_handle,'xlim');
ylim = get(axis_handle,'ylim');
axis_position = get(axis_handle,'position');


%calculate width of figure to be printed
xrange = xlim(2) - xlim(1);
axis_width = axis_position(3);
axis_paperwidth = xrange*(1/x_per_paperinch);
fig_paperwidth = axis_paperwidth*fig_width/axis_width;


%calculate height of figure to be printed
yrange = ylim(2) - ylim(1);
axis_height = axis_position(4);
axis_paperheight = yrange*(1/y_per_paperinch);
fig_paperheight = axis_paperheight*fig_height/axis_height; 


% set paper parameters
set(fig_handle, 'papertype', papertype);
set(fig_handle, 'paperorientation', orient);
%for current paper setup, get size
papersize = get(fig_handle, 'PaperSize');
%set figure position
set(fig_handle, 'PaperPosition', [left_margin, papersize(2)-top_margin-fig_paperheight, fig_paperwidth, fig_paperheight]);

%reset unit values
set(fig_handle,'PaperUnits',t_paperunits);
set(fig_handle,'Units',t_fig_units);
set(axis_handle,'Units',t_axis_units);

