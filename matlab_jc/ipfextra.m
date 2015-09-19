function ipfextra(n,h)

% Tom O'Haver,  Version 1.3, October 23, 2006.

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
global extra


extra=n;

FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start);

axes(h);
h2=gca;