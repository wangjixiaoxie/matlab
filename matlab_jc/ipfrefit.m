function ipfrefit(n,h)
% Computes fit and re-draws graph when the Re-fit slider is moved
% Tom O'Haver, October 2006

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

startnow=start;
delta=n*(max(xx)-min(xx))/100;
for k=1:2:2*NumPeaks,
    startnow(k)=start(k)+randn*delta;
end
FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,startnow);
axes(h);
h2=gca;