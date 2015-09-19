function [pitch_data,ifdgram]=jc_PitchData610avg(data,N,OVERLAP,sigma,F_low,F_high,harmonic)
%weighted avg of several harmonics instead of a single harmonic

%Parameters that one may want to change.
SAMPLING=32000;       
ZOOMT=1;
ZOOMF=3;
TL=5;
FL=5;

[ifdgram,sonogram]=ifdv1(data,SAMPLING,N,OVERLAP,sigma,ZOOMT,ZOOMF,TL,FL); 
%The following is some processing of the output from Tim Gardner's
%program that facilitates indexing and calculations in the rest of the
%program.
a=size(ifdgram);
total_bins=a(2);
Nyquist=SAMPLING/2;
step=Nyquist/a(1); %This is the width (in ms) between bins.
mini=round(F_high/step);
maxi=round(F_low/step);
  
%This loop determines the pitch (fundamental frequency) at each time bin by
%taking the peak of the autocorrelation function within a range specified
%by the user.
for currenttime_bin=1:total_bins;
    slice=ifdgram(:,currenttime_bin); %power at each frequency for the current time bin
    Powerful(1)=0;
    Freqbinest(1)=0;
    for i=1:harmonic %:highest_harmonic
       minimum=a(1)-mini*i;  
       maximum=a(1)-maxi*i;
       freq_window=slice(minimum:maximum);
       [Pow,Ind]=max(freq_window);  
       % Interpolation--subtraction of 0.5 because we care about the
       % central frequency of the bin, not the upper boundary.
       if Ind==1 || Ind==length(freq_window)
          Indest=Ind;
       else
          Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
       end
       Real_Index=minimum+Indest-2; % -1 to account for window; -0.5 (-1?!) to account for the center of the bin
       Freqbinest(i)=(a(1)-Real_Index)/i;
       Powerful(i)=Pow;
    end
    normalizer=sum(Powerful);
    normpower=Powerful/normalizer;
    freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
end
pitch_data=(freqbin_estimate)*step;
