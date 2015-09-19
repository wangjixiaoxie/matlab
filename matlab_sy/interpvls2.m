function [xout,yout]=interpvls2(xvls,yvls)
xout=xvls;
yout=yvls;
ct=1;
if(length(xvls)>1)
    for ii=2:length(xvls)
        xdiff=xvls(ii)-xvls(ii-1);
        ydiff=yvls(ii)-yvls(ii-1);
        xpst=xvls(ii);
        ypst=yvls(ii);
        xpstall=xout(ct+1:end);
        xpstall=makerow(xpstall);
        ypstall=yout(ct+1:end);
        xpre=xvls(ii-1);
        ypre=yvls(ii-1);
        xpreall=xout(1:ct);
        xpreall=makerow(xpreall);
        ypreall=yout(1:ct);
        if(xdiff>1)
           for jj=1:xdiff-1
               addx(jj)=xpre+jj
               addy(jj)=ypre+jj*((ydiff/(xdiff+1)))  
           end
         ct=xpre+xdiff;
         xout=[xpreall addx xpstall]
         yout=[ypreall addy ypstall]
        else
            ct=ct+1;
        end
        
        
    end
else
end