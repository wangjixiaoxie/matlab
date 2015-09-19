function htext=lt_plot_pvalue(p)

%% LT 5/13/15 - plots p value as text
% p = 0.02 (give a number

Ylim=ylim;
Xlim=xlim;

htext=text(Xlim(2)-(Xlim(2)-Xlim(1))/5, Ylim(2)-(Ylim(2)-Ylim(1))/10, ['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
if p<0.05;
    plot(Xlim(2)-(Xlim(2)-Xlim(1))/5, Ylim(2)-2*(Ylim(2)-Ylim(1))/10,'*','MarkerSize',7,'Color','r');
end
