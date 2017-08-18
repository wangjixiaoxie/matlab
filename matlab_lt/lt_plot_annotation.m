%% 11/14/15 - runs matlab annotation, but with my defaults for text box
function lt_plot_annotation(position, text, color);
% e.g. inputs
% position = 1 (top left)
% text = 'test'
% color= 'k'

% if have multiple entries, then put text and color in cell [ COLOR DOES
% NOT WORK YET - picks first color for all text
if ~exist('color','var')
    color='k';
end
    

pos=get(gca, 'Position');


switch position
    case 1
h=annotation('textbox', [pos(1)+0.05 pos(2)+pos(4)-0.15 0.15 0.02], 'String', text);
    case 2
h=annotation('textbox', [pos(1)+0.05 pos(2)+0.15 0.15 0.02], 'String', text);
    case 3
h=annotation('textbox', [pos(1)+0.5 pos(2)+0.15 0.15 0.02], 'String', text);
    case 4
h=annotation('textbox', [pos(1)+0.5 pos(2)+pos(4)-0.15 0.15 0.02], 'String', text);
end


if iscell(color)
set(h, 'FontSize', 13, 'FontWeight', 'bold', 'Color', color{1})
else
  set(h, 'FontSize', 13, 'FontWeight', 'bold', 'Color', color);
end
  

