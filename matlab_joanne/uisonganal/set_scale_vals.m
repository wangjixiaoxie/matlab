function set_scale_vals(fig_handle, x_per_axis_inch, y_fig_inches)

%sets values of the selected figure and its current axis to specified scale for screen display
% x_per_axis_inch = number of units along the x axis per inch of displayed figure
% y_fig_inches = total length of y dimension of displayed figure

% will need the handles for both the spectrogram and the amplitude display
global ih_main_amp ih_main_spect ih_main_labels
userdata = get(fig_handle,'userdata');
h_main_amp = userdata(ih_main_amp);
h_main_spect = userdata(ih_main_spect);
h_main_labels = userdata(ih_main_labels);

%save figure unit values and set to inches
t_fig_units = get(fig_handle,'Units');
set(fig_handle,'Units','inches');

t_axis_units = get(h_main_amp,'Units');
set(h_main_amp,'Units','normalized');  %so that axes rescale appropriately when figure is scaled

%set y axis parameters
fig_position = get(fig_handle,'Position'); 
fig_position(4) = y_fig_inches;
set(fig_handle, 'Position', fig_position);


%set x_axis scaling
set(h_main_amp,'Units','inches');
axis_position = get(h_main_amp,'position');
x_axis_inches = axis_position(3);
x_axis_range = x_per_axis_inch*x_axis_inches;
xlim = get(h_main_amp,'xlim');
xlim(2)=xlim(1)+x_axis_range;
set(h_main_amp,'xlim',xlim);
set(h_main_spect,'xlim',xlim);
set(h_main_labels,'xlim',xlim);

%reset unit values
set(fig_handle,'Units',t_fig_units);
set(h_main_amp,'Units',t_axis_units);
