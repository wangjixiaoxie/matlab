function [time_axis,pitches]=jc_PitchCurve(s,N,OVERLAP,sigma,F_low,F_high)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchCurve(data,1024,1010,2)


%Parameters that one may want to change.
SAMPLING=44100;       
ZOOMT=1;
ZOOMF=3;
TL=5;
FL=5;

%Calls Tim Gardner's program, which takes the phase into
%consideration and perhaps even uses consensus of different sampling rates
%to determine the instantaneous fourier transform.
[ifdgram,sonogram]=ifdv(s,SAMPLING,N,OVERLAP,sigma,ZOOMT,ZOOMF,TL,FL);  

%The following is some processing of the output from Tim Gardner's
%program that facilitates indexing and calculations in the rest of the
%program.
a=size(ifdgram);
total_bins=a(2);
Nyquist=SAMPLING/2;
step=Nyquist/a(1); %This is the width (in ms) between bins.
mini=round(a(1)+F_low/step);
maxi=round(a(1)+F_high/step);

%This loop determines the pitch (fundamental frequency) at each time bin by
%taking the peak of the autocorrelation function within a range specified
%by the user.
for current_bin=1:total_bins
    this_slice=ifdgram(:,current_bin);
    c=xcorr(this_slice);
    freq_window=c(mini:maxi);
    [Fmax,Index]=max(freq_window);
    frequencybins(current_bin)=Index+mini-1-a(1);
end
pitches=frequencybins*step;

% This loop converts the bin number to time for use as an
% x-axis.
for count=1:total_bins
    time_axis(count)=(1/SAMPLING)*(length(s)/total_bins)*count;
end

