function lt_getfig(fig_numbers)

for f=fig_numbers;
    figure_name=['figure(' num2str(f) ')'];
    openfig(figure_name);
end
end

    