%plotshortpulse.m

function []=plotshort(outplot);
%plot downward lines as lines from 3 to 3+std down
%plot upward lines as lines form 3 to 3-std down.
figure
%x value
col={'k' 'r' 'b' 'g'}
lw=[4 2]

for ii=1:4

    xvls=outplot(ii).tms
    yvlinit=(outplot(ii).ct-outplot(ii).bas)/outplot(ii).std
    yvlfin=yvlinit+outplot(ii).off/outplot(ii).std
    ind{1}=find(outplot(ii).st>20)
    ind{2}=find(outplot(ii).st<25)
    er=outplot(ii).er
    makeplots(xvls,yvlinit,yvlfin, er, ind,lw,col{ii},outplot(ii).std);


end

function []=makeplots(xvlin,yvlinit,yvlfin, er, ind, lw,col,sd)
    for ii=1:2
        crind=ind{ii}
        xvls=xvlin(crind);
        yinit=yvlinit(crind);
        yfin=yvlfin(crind);
        crer=er(crind)./sd;
        plot([xvls ;xvls],[yinit ;yfin],'Color',col,'Linewidth',lw(ii))
        hold on;
        plot([xvls-2.5;xvls+2.5],[yfin-crer;yfin-crer],'Color',col,'Linewidth',lw(ii));
        plot([xvls-2.5;xvls+2.5],[yfin+crer;yfin+crer],'Color',col,'Linewidth',lw(ii));
    end
