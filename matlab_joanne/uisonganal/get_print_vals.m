function [orient, papertype, x_per_paperinch, y_per_paperinch] = get_print_vals(fig_handle)

%gets values from the selected figure and its current axis for printing
% orient = paper orientation
% papertype = paper type
% x_per_paperinch = number of units along the x axis per inch of printed figure
% y_per_paperinch = number of units along the y axis per inch of printed figure

%get axis that is curent within the selected figure
axis_handle = get(fig_handle,'CurrentAxes');

%save unit values and set to inches
t_paperunits = get(fig_handle,'PaperUnits');
set(fig_handle,'PaperUnits','inches');
t_fig_units = get(fig_handle,'Units');
set(fig_handle,'Units','inches');
t_axis_units = get(axis_handle,'Units');
set(axis_handle,'Units','inches');

%get other printing parameters
orient = get(fig_handle,'PaperOrientation');
papertype = get(fig_handle,'PaperType');
paperposition = get(fig_handle,'PaperPosition');
fig_position = get(fig_handle,'Position'); 
fig_width = fig_position(3);
fig_paperwidth = paperposition(3);
fig_height = fig_position(4);
fig_paperheight = paperposition(4);

xlim = get(axis_handle,'xlim');
ylim = get(axis_handle,'ylim');
axis_position = get(axis_handle,'position');


%calculate inches per unit of x axis
xrange = xlim(2) - xlim(1);
axis_width = axis_position(3);
axis_paperwidth = (axis_width/fig_width)*fig_paperwidth;
x_per_paperinch = xrange/axis_paperwidth; 
 

%calculate inches per unit of y axis
yrange = ylim(2) - ylim(1);
axis_height = axis_position(4);
axis_paperheight = (axis_height/fig_height)*fig_paperheight;
y_per_paperinch = yrange/axis_paperheight; 


%reset unit values
set(fig_handle,'PaperUnits',t_paperunits);
set(fig_handle,'Units',t_fig_units);
set(axis_handle,'Units',t_axis_units);

