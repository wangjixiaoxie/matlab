function [xout,yout]=interpvls3(xvls,yvls)
xout=xvls;
yout=yvls;
ct=1;
if(length(xvls)>1)
    for ii=2:length(xvls)
        %current margin between consecutive values
        xdiff=xvls(ii)-xvls(ii-1);
        ydiff=yvls(ii)-yvls(ii-1);
        
        xpst=xvls(ii);
        ypst=yvls(ii);
        
        
        xpstall=xvls(ii:end);
        xpstall=makerow(xpstall);
        ypstall=yvls(ii:end);
        
        xpre=xvls(ii-1);
        ypre=yvls(ii-1);
        
        xpreall=xout(1:ct);
        xpreall=makerow(xpreall);
        ypreall=yout(1:ct);
        
        if(xdiff>1)
           for jj=1:(xdiff-1)
               addx(jj)=xpre+jj
               addy(jj)=ypre+jj*((ydiff/(xdiff)))  
           end
         ct=ct+xdiff;
         xout=[xpreall addx xpstall]
         yout=[ypreall addy ypstall]
         clear addx addy
        else
            ct=ct+1;
        end
        
        
    end
else
end