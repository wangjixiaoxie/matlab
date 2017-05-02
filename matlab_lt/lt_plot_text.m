%% LT 8/13/15 - just like text() but with default font size and weight
function lt_plot_text(x, y, string, color, fontsize)

if ~exist('color','var');
    color='k';
end

if ~exist('fontsize','var');
    fontsize=13;
end

text(double(x), double(y), string, 'FontSize', fontsize, 'FontWeight', 'bold', 'Color', color);

