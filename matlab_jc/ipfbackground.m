function ipfbackground(n,h)
% Allows user to click graph to enter background points, 
% then computes fit and re-draws graph
% Tom O'Haver,  Version 1.2, October 17, 2006.

global x
global y
global xx
global yy
global NumPeaks
global Shape
global delta
global c
global FitResults
global start

% Acquire background points from user mouse clicks
subplot(2,1,1);xlabel('Click on the background to the LEFT the peak(s).')
[X1,Y1] = GINPUT(1);
subplot(2,1,1);xlabel('Now click on the background to the RIGHT the peak(s).')
[X2,Y2] = GINPUT(1);
n=length(xx);

% Create "start" for this number of peaks
yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);


FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start);

% axes(h);
h2=gca;