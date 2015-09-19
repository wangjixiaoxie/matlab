function ipfpeaks(n,h)
% Updates "NumPeaks" and re-draws graph when the # Peaks slider is moved
% T. C. O'Haver (toh@umd.edu),     Version 1.6, November 10, 2006.

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

ShapeString=''; 
if round(n)~=NumPeaks,
    NumPeaks=round(n);
    switch Shape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';    
    case 3
        ShapeString='logistic';
    case 4
        ShapeString='Pearson7';
    case 5        
        ShapeString='ExpGaussian';
    otherwise
        ShapeString='';
    end
  subplot(2,1,2)
  xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
end

n=max(xx)-min(xx);
% Create a start value for this number of peaks
start=[];
startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
for marker=1:NumPeaks,
    markx=startpos(marker);
    start=[start markx n/(5*NumPeaks)];
end
RedrawSignal(x,y,xo,dx);
  xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
% FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start);
axes(h);
h2=gca;