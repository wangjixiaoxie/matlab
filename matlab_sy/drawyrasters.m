%draw lines between certain x values at a y height.


function []=drawyrasters(xval,ypts,col)
xvals=[xval;xval];
yval=[ypts(1)*ones(size(xval));ypts(2)*ones(size(xval))];
plot(xvals,yval,col,'Linewidth',1)
