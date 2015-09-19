function [pitch_data]=jc_PitchData729_1024w(shifted,N,OVERLAP,sigma,F_low,F_high,window,shift)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

m=size(shifted);
NumberOfNotes=m(1);

%Parameters that one may want to change.
SAMPLING=44100;      

for found_note=1:NumberOfNotes
    current_note=shifted(found_note,:);
    [sonogram]=sonodvA(current_note,SAMPLING,N,OVERLAP,sigma); 
    %The following is some processing of the output from Tim Gardner's
    %program that facilitates indexing and calculations in the rest of the
    %program.
    a=size(sonogram);
    total_bins=a(2);
    Nyquist=SAMPLING/2;
    step=Nyquist/a(1); %This is the width (in ms) between bins.
    mini=round(F_high/step);
    maxi=round(F_low/step);
    highest_harmonic=2; 
    
    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user.
    for currenttime_bin=1:total_bins
        slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:highest_harmonic %:highest_harmonic %:highest_harmonic
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
            Real_Index=minimum+Indest-1; % -1 to account for window;
            Freqbinest(i)=(a(1)-Real_Index)/i;
            Powerful(i)=Pow;
        end
        normalizer=sum(Powerful);
        normpower=Powerful/normalizer;
        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
    end
    pitch_data(found_note).pitches=freqbin_estimate*(Nyquist/(a(1)-1));
end
        
%Plotting - like pplotter
figure; hold on;
for i=1:15
    plotshiftedD=pitch_data(i).pitches+180*(i);
    plot(plotshiftedD)
end
        
        