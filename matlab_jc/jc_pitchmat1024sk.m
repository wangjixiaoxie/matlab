function [pitch_data,avg]=jc_pitchmat1024sk(shifted,N,OVERLAP,sigma,F_low,F_high,harms,filetype,returnAV)

% shifted - each row is the raw data for one rendition of the syllable of
    % interest - note that you lose the first and last 1024points
% N - default is 1024
% Overlap - default is 1020
% sigma - based on the spectral sparseness - so it should be calibrated to
    % the fundamental frequency b/c that is the distance between harmonics -
    % for ZF song (with harmonics ~ 750Hz), I use about sigma=2
% filetype - 'w' if wav files


%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%jc_pitchmat1024(shifted,1024,1010,2,2000,2700,[1 2 3],'obs0');

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

%return the average note
if returnAV==1
    for i=1:size(pitch_data,1);avg(i)=mean(pitch_data(i,:));end
end
% Note that this method of calculating the average agrees with the previous
% method over the smooth part of the note but may be useful for eliminating
% transients.
%         HI=F_high-20;
%         LO=F_low+20;
% if returnAV==1
%     for i=1:size(pitch_data,1)
%         ab=pitch_data(i,:);
%         summer=0;
%         count=0;
%         for j=1:length(ab)
%             if F_low+30<ab(j)&& ab(j)<F_high-30
%                 summer=summer+ab(j);
%                 count=count+1;
%             end
%         end
%         if count==0
%             avg(i)=mean(ab);
%         else
%         avg(i)=summer/count;
%         end
%     end
% end
