function ipfzoom(n,h)
% Updates "start" and re-draws graph when the Zoom slider is moved
% Tom O'Haver, October 2006

global x
global y
global xx
global yy
global xo
global dx
global start

dx=round(n);
[xx,yy,start]=RedrawSignal(x,y,xo,dx);

axes(h);
h2=gca;