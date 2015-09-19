function ipfpeaksmouseup
% Computes fit on a mouseup in the # Peaks slider
% T. C. O'Haver (toh@umd.edu),   test  Version 1.6, November 10, 2006.

global x
global y
global xx
global yy
global xo
global dx
global NumPeaks
global Shape
global delta
global c
global FitResults
global start

  
FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start);
% axes(h);
h2=gca;