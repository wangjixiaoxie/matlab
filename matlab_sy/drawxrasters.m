%draw lines between certain x values at a y height.

function []=drawxrasters(xval,xval2,yarray,col)

xvals=[xval;xval2];
yval=[yarray*ones(size(xval));yarray*ones(size(xval))];
plot(xvals,yval,col,'Linewidth',3)
 