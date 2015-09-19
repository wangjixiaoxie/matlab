function [pitch_data]=jc_PitchDataVCOifd(x,N,OVERLAP,sigma,F_low,F_high)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

%Parameters that one may want to change.
SAMPLING=32000;      

    %[sonogram]=sonodv(x,SAMPLING,N,OVERLAP,sigma); 
    [sonogram]=ifdv(x,SAMPLING,N,OVERLAP,sigma,1,15,5,5);
    %The following is some processing of the output from Tim Gardner's
    %program that facilitates indexing and calculations in the rest of the
    %program.
    a=size(sonogram);
    total_bins=a(2);
    Nyquist=SAMPLING/2;
    step=Nyquist/a(1); %This is the width (in ms) between bins.
    mini=round(F_high/step);
    maxi=round(F_low/step);
    
    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user.
    for currenttime_bin=1:total_bins
        slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:1 %:highest_harmonic %:highest_harmonic
            minimum=a(1)-mini*i;  
            maximum=a(1)-maxi*i;
            freq_window=slice(minimum:maximum);
            [Pow,Ind]=max(freq_window);  

            Real_Index=minimum+Ind-1.25; % -1 to account for window;(-1?!) to account for weird error in ifdv
            Freqbinest(i)=(a(1)-Real_Index)/i;
            Powerful(i)=Pow;
        end
        normalizer=sum(Powerful);
        normpower=Powerful/normalizer;
        freqbin_estimate(currenttime_bin)=Freqbinest(i);%dot(normpower,Freqbinest);
    end
    pitch_data=freqbin_estimate*step;
        
        
        
        