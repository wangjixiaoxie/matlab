function assignlabel(varargin)
% ax = get(ploth,'parent');
% set(ploth,'color','g')
% ud = get(ploth,'UserData');
fh = varargin{1};
key = varargin{2}.Key;
ud = get(fh,'UserData');
ud{1} = key;
set(fh,'UserData',ud)

ax = get(fh,'Children');
set(ax(end),'color','g')
axes(ax(end))
text(0.3, 0.7, key, 'FontSize',15)
% lh = ud{1};
% % flagformoving = ud{2};
% flagformoving(lh) = 1;
% ud{2} = flagformoving;
% set(ax,'UserData',ud)
end