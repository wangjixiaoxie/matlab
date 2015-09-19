function [pitch_data,ifd,sono]=jc_PitchData606p(N,OVERLAP,sigma,F_low,F_high,data)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

%Sets parameters and runs findwnote4 to extract all instances of the
%desired note in the specified context.
Timeshift=0;
Note_length=6500;
UseFit=1;
ChanSpec='obs0';
%[fvalsstr]=findwnote4(batchFWN,noteFWN,prenoteFWN,postnoteFWN,Timeshift,[F_low F_high],Note_length,UseFit,ChanSpec);

%Parameters that one may want to change.
SAMPLING=32000;       
ZOOMT=1;
ZOOMF=3;
TL=5;
FL=5;


for found_note=1:1 %length(fvalsstr)
    
    %Calls Tim Gardner's program, which takes the phase into
    %consideration and perhaps even uses consensus of different sampling rates
    %to determine the instantaneous fourier transform.
    %current_note=fvalsstr(found_note).datt;
    %chopped_note=current_note(1000:5500);
    [ifdgram,sonogram]=ifdv(data,SAMPLING,N,OVERLAP,sigma,ZOOMT,ZOOMF,TL,FL); 
    ifd(found_note).gram=ifdgram;
    sono(found_note).gram=sonogram;
    %y(found_note)=sonogram;
    %The following is some processing of the output from Tim Gardner's
    %program that facilitates indexing and calculations in the rest of the
    %program.
    a=size(ifdgram);
    total_bins=a(2);
    Nyquist=SAMPLING/2;
    step=Nyquist/a(1); %This is the width (in ms) between bins.
    mini=round(F_high/step);
    maxi=round(F_low/step);
    highest_harmonic=round(a(1)/mini)-1;  %This often throws away the last one, but who cares
    
    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user.
    for currenttime_bin=1:total_bins
        slice=ifdgram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:highest_harmonic
            minimum=a(1)-mini*i+1;  % Add one so that the freq_window widths at different harmonics are integer multiples
            maximum=a(1)-maxi*i;
            freq_window=slice(minimum:maximum);
            [Pow,Ind]=max(freq_window);  
            % Interpolation
            if Ind==1 || Ind==length(freq_window)
                Indest=Ind;
            else
                Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
            end
            Freqbinest(i)=mini-(Indest-0.5)/i; %Subtracting 0.5 accounts for the bottom and top of the slice.
            Powerful(i)=Pow;
        end
        normalizer=sum(Powerful);
        normpower=Powerful/normalizer;
        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
    end
    pitch_data(found_note).pitches=freqbin_estimate*step;
end

        
        
        
        
        
        