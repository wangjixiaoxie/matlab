function [h_print] = print_song

%calls uiprint to print the current figure window at current figure dimensions

%global variables for margins
global pr_left_margin pr_bottom_margin 


fig_handle = gcf;

%save unit values and set to inches
t_paperunits = get(fig_handle,'PaperUnits');
set(fig_handle,'PaperUnits','inches');
t_fig_units = get(fig_handle,'Units');
set(fig_handle,'Units','inches');

position = get(fig_handle,'position');
position(1) = pr_left_margin;
position(2) = pr_bottom_margin;
set(fig_handle,'paperposition',position);

%reset unit values
set(fig_handle,'PaperUnits',t_paperunits);
set(fig_handle,'Units',t_fig_units);



%call print dialog box
%uiprint_sa;
printdlg;
