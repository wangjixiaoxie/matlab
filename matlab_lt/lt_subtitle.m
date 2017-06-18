%% LT 4/29/15 - like subtitle, but fonts formatted
function lt_subtitle(string)

% [~, h]=subtitle(string);
% set(h,'FontSize',16,'FontWeight','bold')

h=mtit(string);

pos=get(h.th, 'position');
set(h.th,'FontSize',16,'FontWeight','bold','position',[pos(1) pos(2)+0.05 pos(3)])