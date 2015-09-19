function ipfzoomup
% ZOOMS UP one point and redraws graph
% Tom O'Haver,   Version 1.6, November 10, 2006.

global x
global y
global xx
global yy
global xo
global dx
global start

dx=dx+1;
[xx,yy,start]=RedrawSignal(x,y,xo,dx);

% axes(h);
h2=gca;