function ipfcustom(n,h)
% Allows user to click graph to enter start positons, 
% then computes fit and re-draws graph
% Tom O'Haver,    Version 1.6, November 10, 2006.

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

% Acquire first-guess peak positions from user mouse pointer
subplot(2,1,1);xlabel('Click on the estimated positions of each proposed component peak.')
[clickX,clickY] = GINPUT(NumPeaks);

% Create a "start" vector using these peak positions, with peak widths equal 
% to 1/10 of the zoom region.
n=max(xx)-min(xx);
width=n/(5*NumPeaks);
start=[];
for k=1:NumPeaks,
    start=[start clickX(k) width];
end

FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start);

% axes(h);
h2=gca;