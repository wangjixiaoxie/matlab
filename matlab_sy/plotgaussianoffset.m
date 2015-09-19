function plotgaussianoffset(axh,xbnds,mu,sigma,color,lstyle,flip)
axes(axh);
pt1=1/(sigma*sqrt(2*pi))
pt2a=exp(-((xbnds-mu(1)).^2)/(2*(sigma^2)));
pt2b=exp(-((xbnds-mu(2)).^2)/(2*(sigma^2)));
y1=pt1*pt2a*[xbnds(2)-xbnds(1)];
y2=pt1*pt2b*[xbnds(2)-xbnds(1)];

y=(y1+y2)/2;


if(flip)
   if(exist('lstyle'))
        plot(y,xbnds,'Color',color,'Linestyle',lstyle,'Linewidth',4);
    else
        plot(y,xbnds,'Color',color);
    end    
else
    if(exist('lstyle'))
        plot(xbnds,y,'Color',color,'Linestyle',lstyle,'Linewidth',4);
    else
        plot(xbnds,y,'Color',color);
    end
end