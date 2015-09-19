function [xx,yy,start]=RedrawSignal(x,y,xo,dx);
% Plots the entire signal (x,y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom sliders.
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008.
global NumPeaks
% Top half of the figure shows original signal
Startx=round(xo-(dx/2));
Endx=round(xo+(dx/2)-1);
if Endx>length(y),
    Endx=length(y);
end
if Startx<1,
     Startx=1;
end
PlotRange=[Startx:Endx];
if PlotRange<5, PlotRange=[xo:xo+5];,end
xx=x(PlotRange);
yy=y(PlotRange); 

% Comment out the next 5 lines to disable auto-zero zoom operation
% X1=min(xx);
% X2=max(xx);
% Y1=yy(1);
% Y2=yy(length(yy));
% yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);

subplot(2,1,1);plot(xx,yy,'b.'); % Plot the original signal in blue
hold on

% Mark starting peak positions with vertical dashed lines
% Determine locations for peak (vertical line) markers
n=max(xx)-min(xx);
start=[];
startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
for marker=1:NumPeaks,
    markx=startpos(marker);
    start=[start markx n/5];
    subplot(2,1,1);plot([markx markx],[min(yy) max(yy)],'m--')
end
hold off
title('Use Pan and Zoom sliders to isolate peaks to be fit in upper window.')
xlabel('Line up the estimated peak positions roughly with the vertical lines')
% axis([round(xo-(dx/2)) round(xo+(dx/2)-1) min(yy) max(yy)])
axis([x(Startx) x(Endx) min(yy) max(yy)])

% Bottom half of the figure shows full signal
subplot(2,1,2);plot(x,y)
hold on
for marker=1:NumPeaks,
    markx=startpos(marker);
    subplot(2,1,2);plot([markx markx],[min(y) max(y)],'m--')
end
hold off
xlabel(['To fit peaks in upper window, set # Peaks and Shape, then click Re-fit.'])