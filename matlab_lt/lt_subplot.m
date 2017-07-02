function h=lt_subplot(arg1, arg2, arg3);
%% LT 5/4/15 - subplot, formatted as I want

h=subplot(arg1, arg2, arg3);

lt_plot_format;

% subplotsqueeze(gca, 1.2);
end
