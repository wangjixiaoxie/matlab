function ipfnudgedown
% Pans down one point and redraws graph
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008.

global x
global y
global xx
global yy
global xo
global dx
global start

xo=xo-1;
if xo<1,xo=1,end
[xx,yy,start]=RedrawSignal(x,y,xo,dx);
% axes(h);
h2=gca;