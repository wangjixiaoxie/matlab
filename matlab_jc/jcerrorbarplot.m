function jcerrorbarplot(x,y,semx,semy,clr)
plot(x,y,'Marker','.','Color',clr,'Markersize',15)
plot([x x],[y-semy y+semy],'Color',clr)
plot([x-semx x+semx],[y y],'Color',clr)