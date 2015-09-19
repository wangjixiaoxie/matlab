function ipfshape(n,h)
% Redraws signal when the shape slider is moved
% T. C. O'Haver (toh@umd.edu),    Version 1.6, November 10, 2006.

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

if round(n)~=Shape,
    Shape=round(n);
        switch Shape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';    
    case 3
        ShapeString='logistic';
    case 4
        ShapeString='Pearson';
    case 5        
        ShapeString='ExpGaussian';
    otherwise
    end
  subplot(2,1,2)
  xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
end
% FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start);
axes(h);
h2=gca;