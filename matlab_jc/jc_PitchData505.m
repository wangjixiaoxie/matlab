function [pitch_data,ifd,sono]=jc_PitchData505(N,OVERLAP,sigma,F_low,F_high,batchFWN,noteFWN,prenoteFWN,postnoteFWN)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

%Sets parameters and runs findwnote4 to extract all instances of the
%desired note in the specified context.
Timeshift=0;
Note_length=4096;
UseFit=1;
ChanSpec='obs0';
[fvalsstr]=findwnote4(batchFWN,noteFWN,prenoteFWN,postnoteFWN,Timeshift,[F_low F_high],Note_length,UseFit,ChanSpec);

%Parameters that one may want to change.
SAMPLING=32000;       
ZOOMT=1;
ZOOMF=3;
TL=5;
FL=5;


for found_note=1:length(fvalsstr)
    
    %Calls Tim Gardner's program, which takes the phase into
    %consideration and perhaps even uses consensus of different sampling rates
    %to determine the instantaneous fourier transform.
    current_note=fvalsstr(found_note).datt;
    [ifdgram,sonogram]=ifdv(current_note,SAMPLING,N,OVERLAP,sigma,ZOOMT,ZOOMF,TL,FL); 
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
    mini=round(a(1)+F_low/step);
    maxi=round(a(1)+F_high/step);

    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user.
    for current_bin=1:total_bins
        this_slice=ifdgram(:,current_bin);
        c=xcov(this_slice);
        freq_window=c(mini:maxi);
        [Fmax,Index]=max(freq_window);
        if Index==1 || Index==length(freq_window) 
            peak=Index;
        else
            peak=pinterp([Index-1;Index;Index+1],[freq_window(Index-1);freq_window(Index);freq_window(Index+1)]);
        end
        frequencybins(current_bin)=peak+mini-1-a(1);
    end
    pitch_data(found_note).pitches=frequencybins*step;
end

    