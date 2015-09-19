function lt_global_subplot_title(title_string);
% LT 9/11/14 - for subplot, will add single title on top of everything.

[ax, h]=subtitle(title_string);
set(h,'Position',title_pos);

end