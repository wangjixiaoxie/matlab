function ipfpan(n,h)
% Re-draws graph and updates "start" when the Pan slider is moved
% Tom O'Haver, July 2006

global x
global y
global xx
global yy
global xo
global dx
global start

xo=round(n);
[xx,yy,start]=RedrawSignal(x,y,xo,dx);
axes(h);
h2=gca;