function htext=lt_plot_pvalue(p,string,position)
% position = 1, 2 (corner of plot)

if ~exist('position','var');
    position = 1;
end

    
%% LT 5/13/15 - plots p value as text
% p = 0.02 (give a number

Ylim=ylim;
Xlim=xlim;

if position ==1;
xpos=Xlim(2)-(Xlim(2)-Xlim(1))/5;
elseif position==2;
xpos=Xlim(1)+(Xlim(2)-Xlim(1))/5;
end

htext=text(xpos, Ylim(2)-(Ylim(2)-Ylim(1))/10, ['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
if p<0.05;
    plot(xpos, Ylim(2)-2*(Ylim(2)-Ylim(1))/10,'*','MarkerSize',7,'Color','r');
end

try
    lt_plot_text(xpos, Ylim(2)-3*(Ylim(2)-Ylim(1))/10, string);
end
