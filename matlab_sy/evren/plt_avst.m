function plt_avst(tt,av,st,clr);
% plt_avst(tt,av,st,clr);
%
fill([tt,tt(end:-1:1)],[av+st,av(end:-1:1)-st(end:-1:1)],clr,'EdgeC',clr);
return;
