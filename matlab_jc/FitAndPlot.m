function FitResults=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start)
% Given isolated segment of signal (xx,yy), plots it in upper half, computes fit with
% "NumPeaks" component peaks of shape "Shape", starting with start values
% "start", then plots residuals in lower half. 
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008, uses unconstrained fit.

global c
global extra

n=length(xx);
x0=min(xx);

% Set lower bound to zero (avoids negative peak positions or widths)
LB=zeros(size(1:(2*NumPeaks)));
UB=[]; % No upper bound on peak positions or widths
h=figure(1);

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.01);
switch Shape
    case 1
        FitParameters=FMINSEARCH('fitgaussian',start,options,xx,yy);
        ShapeString='Gaussian';
    case 2
        FitParameters=FMINSEARCH('fitlorentzian',start,options,xx,yy);
        ShapeString='Lorentzian';    
    case 3
        FitParameters=FMINSEARCH('fitlogistic',start,options,xx,yy);
        ShapeString='Logistic';
    case 4
        FitParameters=FMINSEARCH('fitpearson',start,options,xx,yy,extra);
        ShapeString='Pearson';
    case 5
        FitParameters=FMINSEARCH('fitexpgaussian',start,options,xx,yy,-extra);
        ShapeString='ExpGaussian';
    otherwise
end

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
       FitResults=[round(m) FitParameters(2*m-1) c(m) FitParameters(2*m)];
    else
       FitResults=[FitResults ; [round(m) FitParameters(2*m-1) c(m) FitParameters(2*m)]];
   end
end

% Construct model from fitted parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
   switch Shape
    case 1
        A(m,:)=gaussian(xx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        A(m,:)=lorentzian(xx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        A(m,:)=logistic(xx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        A(m,:)=Pearson(xx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        A(m,:)=expgaussian(xx,FitParameters(2*m-1),FitParameters(2*m),-extra)';
    otherwise
  end
end
model=c'*A;  % Multiplies each row by the corresponding amplitude and adds them up

% Top half of the figure shows original signal and the fitted model.
subplot(2,1,1);plot(xx,yy,'b.'); % Plot the original signal in blue
hold on
for m=1:NumPeaks,
    plot(xx,c(m)*A(m,:),'g')  % Plot the individual component peaks in green
end
% Mark starting peak positions with vertical dashed lines
for marker=1:NumPeaks,
    markx=start((2*marker)-1);
    subplot(2,1,1);plot([markx markx],[0 max(yy)],'m--')
end
plot(xx,model,'r');  % Plot the total model (sum of compenent peaks) in red
hold off;
axis([min(xx) max(xx) min(yy) max(yy)]);
title('Pan and Zoom to select peaks; # peaks, Shape, Re-fit to fit peaks.')
xlabel('Vertical dotted lines indicate first guess peak positions');

% Bottom half of the figure shows the residuals and displays RMS error between original signal and model
residual=yy-model;
subplot(2,1,2);plot(xx,residual,'b.')
MeanFitError=100*norm(residual)./(sqrt(n)*max(yy));
xlabel(['Extra = ' num2str(extra) '     Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(MeanFitError) '%'] )
axis([min(xx) max(xx) min(residual) max(residual)]);

% Un-comment the next 3 lines to display Fit Results on upper graph
subplot(2,1,1);
residualaxis=axis;
text(residualaxis(1), 0.8*residualaxis(4) ,num2str(FitResults));

% Add "FitResults" and/or "MeanFitError" if you wish to have these printed out
% in the Matlab command window each time a fit is computed.

% FitResults
% MeanFitError