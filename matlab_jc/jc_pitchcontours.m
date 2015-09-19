function [pitch_data]=jc_pitchcontours(fvals,N,OVERLAP,sigma,F_low,F_high,harms,filetype)
%F_low and F_high are the boundary frequencies in which to look for a peak
    %in the autocorrelation function (i.e. in which to calculate the pitch at
    %each point in time).

% sigma of 1.0 is generally good for BF song (with harmonic spacing >>1500Hz)
% sigma of 1.5-2.0 is generally good for ZF song (with harmonic spacing <<1500Hz)
    % this is based on the consensus function from Tim Gardner

% fvals is the output of one of the findwnote functions.  If you don't use
% these, note that each column of shifted (see line ) is about 100ms of
% song aligned to the beginning of the syllable of interest.

% harms allows you to weight the  pitch estimate at each time point based 
% on the power at each harmonic.  Usually, I just choose the 
% single harmonic that tends to have the most power. Either way gives a similar solution.
% Looking only at the first harmonic
    % pitchcontours_pre=jc_pitchcontours(fvals_pre,1024,1020,1,2000,3000,[1],'obs0');
% Looking at all three harmonics below 10000
    % pitchcontours_pre=jc_pitchcontours(fvals_pre,1024,1020,1,2000,3000,[1 2 3],'obs0');
    
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
m=size(shifted);
NumberOfNotes=m(1);

%Parameters that one may want to change.
if strcmp(filetype,'obs0')
    SAMPLING=32000; 
else SAMPLING=44100;
end

%sonodvA program
t=-N/2+1:N/2;
sigma=(sigma/1000)*SAMPLING;
%Gaussian and first derivative as windows.
w=exp(-(t/sigma).^2);
timebins=floor((size(shifted,2)-(N))/(N-OVERLAP))+1;
freqbins=(N/2)+1;

Nyquist=SAMPLING/2;
step=Nyquist/freqbins; %This is the width (in ms) between bins.
mini=round(F_high/step);
maxi=round(F_low/step);
highest_harmonic=1; 

sonogram=zeros(freqbins,timebins);
pitch_data=zeros(timebins,NumberOfNotes);
for found_note=1:NumberOfNotes
    
    sonogram=abs(flipdim(spectrogram(shifted(found_note,:),w,OVERLAP,N,SAMPLING),1)); %gaussian windowed spectrogram

    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user
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
            if Ind==1 || Ind==length(freq_window)
                Indest=Ind;
            else
                Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
            end
            Real_Index=minimum+Indest-1; % -1 to account for window;
            Freqbinest(i)=(freqbins-Real_Index)/harms(i);
            Powerful(i)=Pow;
        end
        normalizer=sum(Powerful);
        normpower=Powerful/normalizer;
        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
    end
    pitch_data(:,found_note)=freqbin_estimate*(Nyquist/(freqbins-1));
end

