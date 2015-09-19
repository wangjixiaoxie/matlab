%% LT 8/13/15 - just like text() but with default font size and weight
function lt_plot_text(x, y, string);

text(x, y,string, 'FontSize', 13, 'FontWeight', 'bold');
