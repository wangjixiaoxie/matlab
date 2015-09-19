% Interactive Peak Fitter for signal in global vectors x,y.
% See http://www.wam.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008, uses unconstrained fit.

global x
global y
global xo
global dx
global NumPeaks
global Shape
global Scan
global Width
global c
global FitResults
global start
global extra

close
format compact
warning off all
%warning('OFF', 'Exiting:')

% Adjust x and y vector shape to 1 x n (rather than n x 1)
x=reshape(x,1,length(x));
y=reshape(y,1,length(y));

% x=[1:length(y)]; % Use this only to create an x vector if needed

%'Type FitResults to see result of the peak fitting'

% Initial values of filter parameters
NumPeaks=1; % Initial Number of peaks in model
Shape=1; % Initial Shape of the peaks (1=Gaussian, 2=Lorentzian, etc)
extra=.5;
delta=0; % Initial Random change in initial start guesses for each re-fit
xo=length(y)/2; % Initial Pan setting
dx=length(y)/4; % Initial Zoom setting

% Plot the signal and its fitted components and residuals
RedrawSignal(x,y,xo,dx);
h=figure(1);
h2=gca;

% Maximum/Minimum ranges of the sliders (change as needed)
MaxPeaks=5; % You can set a higher limit here if you wish (no inherent limit)
MaxShapes=5; % See http://www.wam.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
MaxExtra=10;
MinExtra=.5;
MaxDelta=2;
MaxPan=length(y);
MaxZoom=length(y)/1.1;
MinZoom=length(y)/40;

% Draw the sliders
rtslid(h,@ipfpan,h2,1,'Scale',[1 MaxPan],'Def',xo,'Back',[0.9 0.9 0.9],'Label','Pan','Position',[0.03 0.45 0.03 0.5]);
rtslid(h,@ipfzoom,h2,1,'Scale',[MinZoom MaxZoom],'Def',dx,'Back',[0.9 0.9 0.9],'Label','Zoom','Position',[0.94 0.45 0.03 0.5]);

rtslid(h,@ipfpeaks,h2,0,'Scale',[1 MaxPeaks],'Def',NumPeaks,'Back',[0.9 0.9 0.9],'Label','# Peaks','Position',[0.03 0.30 0.03 0.06]);
rtslid(h,@ipfshape,h2,0,'Scale',[1 MaxShapes],'Def',1,'Back',[0.9 0.9 0.9],'Label','Shape','Position',[0.03 0.17 0.03 0.06]);
rtslid(h,@ipfextra,h2,0,'Scale',[MinExtra MaxExtra],'Def',extra,'Back',[0.9 0.9 0.9],'Label','Extra','Position',[0.03 0.01 0.03 0.09]);

rtslid(h,@ipfrefit,h2,0,'Scale',[0 MaxDelta],'Def',0,'Back',[0.9 0.9 0.9],'Label','Re-fit','Position',[0.94 0.25 0.03 0.12]);
rtslid(h,@ipfcustom,h2,0,'Scale',[0 0],'Def',0,'Back',[0.9 0.9 0.9],'Label','Custom','Position',[0.94 0.14 0.03 0.04]);
rtslid(h,@ipfbackground,h2,0,'Scale',[0 0],'Def',0,'Back',[0.9 0.9 0.9],'Label','BG','Position',[0.94 0.03 0.03 0.04]);

% Attaches KeyPress test function to the figure (for keyboard shortcuts)
KeyPressTest