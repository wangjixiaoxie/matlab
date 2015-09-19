%% LT 7/31/14 - input number of subplots and get out 1) numfigures, 2) rows and 3) columns for a subplot;

function [num_figures, row_plots, col_plots]=lt_get_subplot_size(num_plots,PlotsPerFig)
% INPUTS
% num_plots = 6 (i.e. 6 subplots);
% PlotsPerFig = 9 (i.e. at most 9 plots for each figure);
% PlotsPerFig = 0 (default - will set to 1 figure holding all plots);
% OUTPUTS
% number of figures, and rows/columns in each figure;



if PlotsPerFig==0;
    num_figures=1;
else
    num_figures=ceil((num_plots)/PlotsPerFig);
end


for i=1:num_figures;
    if i==num_figures;
        row_plots(i)=ceil((num_plots-(i-1)*9)/3);
        col_plots(i)=ceil((num_plots-(i-1)*9)/row_plots(i));
    else
        row_plots(i)=3;
        col_plots(i)=3;
    end
end
