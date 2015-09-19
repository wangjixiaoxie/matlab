function [x_per_axis_inch, y_fig_inches] = get_scale_vals(fig_handle)

%gets values from the selected figure and its current axis for scaling screen display
% x_per_axis_inch = number of units along the x axis per inch of displayed figure
% y_fig_inches = total length of y dimension of displayed figure

%get axis that is curent within the selected figure
axis_handle = get(fig_handle,'CurrentAxes');

%save unit values and set to inches
t_fig_units = get(fig_handle,'Units');
set(fig_handle,'Units','inches');
t_axis_units = get(axis_handle,'Units');
set(axis_handle,'Units','inches');

%get parameters
fig_position = get(fig_handle,'Position'); 
y_fig_inches = fig_position(4);

xlim = get(axis_handle,'xlim');
axis_position = get(axis_handle,'position');
x_axis_inches = axis_position(3);
x_per_axis_inch = (xlim(2)-xlim(1))/x_axis_inches;

%reset unit values
set(fig_handle,'Units',t_fig_units);
set(axis_handle,'Units',t_axis_units);

