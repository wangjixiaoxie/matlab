function ipfnudgeup
% Pans up one point and redraws graph
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008.

global x
global y
global xx
global yy
global xo
global dx
global start

ly=length(y);
xo=xo+1;
if xo>ly,xo=ly,end
[xx,yy,start]=RedrawSignal(x,y,xo,dx);
% axes(h);
h2=gca;