% LT 4/25/14 - takes a number (how many categories to plot) and makes that
% many distinct colors, or take one color and make it graded.
% e.g. plot_colors = lt_make_plot_colors(5,0) gives 5 different colors
% e.g. plot_colors = lt_make_plot_colors(5,1,[1 0 0]) gives color [1 0 0],
% but graded equidistantly from 0.25*[1 0 0]
% plot_color is 5-cell cell aray

function plot_color = lt_make_plot_colors(num_col,graded,graded_color)

if graded==1;
    StartCol=0.25*graded_color;
    ColIncrements=(graded_color-StartCol)./num_col;
    for i=1:num_col;
        plot_color{i}=StartCol+ColIncrements.*(i-1);
    end
else
    if num_col==4; % ad hoc, these colors work decent.
        plot_color{1}=[0.1 0.1 0.1];
        plot_color{2}=[0.9 0.1 0.1];
        plot_color{3}=[0.1 0.9 0.1];
        plot_color{4}=[0.1 0.1 0.9];
    else
        split=ceil(num_col/2);
        for i=1:split;
            plot_color{i}=[i/split 0 1] - [0 0 i/split];
        end
        for i=split+1:num_col;
            plot_color{i}=[0.6 (i-split)/(num_col-split+1) 0.5] - [0.6*(i-split)/(num_col-split+1) 0 0.5*(i-split)/(num_col-split+1)];
        end
    end
end