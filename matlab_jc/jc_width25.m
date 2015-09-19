function [sizer]=jc_width25(norm)
x=[0 0 0 125 0 250];
for i=1:251
    xdata(i)=i;
end
sizes=0;
for i=1:20
    ydata=norm(500:750,i);
    ydata=ydata';
    [xf]=lsqcurvefit(@jc_SETgauss,x,xdata,ydata);
    sizes(i)=0.5*xf(1)^2+xf(3)^2+0.5*xf(5)^2;
end
sizer=median(sizes);
    