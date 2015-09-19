function [pitch_data]=jc_PitchData208(shifted,N,OVERLAP,sigma,F_low,F_high,harms,filetype)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

%Sets parameters and runs findwnote4 to extract all instances of the
%desired note in the specified context.

%Parameters that one may want to change.
if strcmp(filetype,'obs0')
    SAMPLING=32000; 
else SAMPLING=44100;
end


for found_note=1:size(shifted,1)
    %Calls Tim Gardner's program, which takes the phase into
    %consideration and perhaps even uses consensus of different sampling rates
    %to determine the instantaneous fourier transform.
    current_note=shifted(found_note,:);
    
    % The value of '30' makes it smoother.
    [sonogram]=ifdvNS(current_note,SAMPLING,N,OVERLAP,sigma,1,5);%,5,5); 

    %The following is some processing of the output from Tim Gardner's
    %program that facilitates indexing and calculations in the rest of the
    %program.
    a=size(sonogram);
%     total_bins=a(2);
    Nyquist=SAMPLING/2;
    step=Nyquist/a(1); %This is the width (in ms) between bins.
    mini=round(F_high/step);
    maxi=round(F_low/step);
%     highest_harmonic=2; %round(a(1)/mini)-1;  %This often throws away the last one, but who cares

    freqbins=a(1);
    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user.
    for currenttime_bin=1:size(sonogram,2)
        slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:length(harms)
            minimum=freqbins-mini*harms(i);  
            maximum=freqbins-maxi*harms(i);
            freq_window=slice(minimum:maximum);
            [Pow,Ind]=max(freq_window);  
            % Interpolation--subtraction of 0.5 because we care about the
            % central frequency of the bin, not the upper boundary.
            if Ind==1 || Ind==length(freq_window) || Ind==0
                Indest=Ind;
            else
                Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
            end
      
            Real_Index=minimum+Indest-1; % -1 to account for window;
            Freqbinest(i)=(freqbins-Real_Index)/harms(i);
            Powerful(i)=Pow;
        end
        normalizer=sum(Powerful);
        if normalizer==0;
            % This is if there is no power
            freqbin_estimate(currenttime_bin)=Freqbinest(1);
        else
            normpower=Powerful/normalizer;
            freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
        end
    end
    pd=freqbin_estimate*(Nyquist/(freqbins-1));
    for i=2:length(pd)-1
        if isnan(pd(i))
            pd(i)=pinterp([i-1;i;i+1], [pd(i-1);pd(i);pd(i+1)]);
        end
    end
    pitch_data(:,found_note)=pd;
end

        
        
        
        
        
        