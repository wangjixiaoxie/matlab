function pitch_data=jc_pitchcontourFV(fvals,N,OVERLAP,sigma,F_low,F_high,harms,filetype)
% LT notes:
% N - size of windows (determines frequency resolution)
% OVERLAP - overlap of successive windows for spectrogram function
% sigma - variance for gaussian window used in the spectrogram function.
    

%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).

% F_low and F_high set the FF range to look for peaks within
% harms are the harmonics to take the weighted average of pitch within
    % (weighted by power at that harmonic)
% sigma of ~1 is optimal for Bengalese song, ~2 is optimal for Zebra song ---
    % the reason for this is the optimal sigma (as determined by the
    % Consensus function in the Gardner/Magnasco papers) is inversely
    % related to the sparseness of the signal (i.e. the width in Hz between
    % harmonics)  --- If your curves for Bengalese song look choppy at
    % sigma=1, the song is probably noisy due to probes/etc. and you should
    % increase the sigma.  There is a rigorous but slow way to do this - 
    % let me know if you want to do it. A decent alternative is just to
    % increase it until the curves look smooth.
    
%jc_pitchcontourFV(fvals,1024,1020,1,2000,2700,[1 2],'obs0');

for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt; % this is sound data
end


NumberOfNotes=size(shifted,1);

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
timebins=floor((size(shifted,2)-(N))/(N-OVERLAP))+1;  % number of bins, i.e. how many times the window N will be shifted to reach the end of the data.
freqbins=(N/2)+1;

Nyquist=SAMPLING/2;
step=Nyquist/freqbins; % This is the width of freqency bins (in Hz)
mini=round(F_high/step);
maxi=round(F_low/step); % step numbers for the min and max of the window you are interested in.

sonogram=zeros(freqbins,timebins);
pitch_data=zeros(timebins,NumberOfNotes);
for found_note=1:NumberOfNotes
    [S,~,~]=spectrogram(shifted(found_note,:),w,OVERLAP,N,SAMPLING);
    sonogram=abs(flipdim(S,1)); % gaussian windowed spectrogram, for entire data duration

    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user
    for currenttime_bin=1:size(sonogram,2)
        slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:length(harms)
            minimum=freqbins-mini*harms(i);  % to extract freq window. harms = 2 just doubles the min and max limits. so the window itself is actually double sized.
            maximum=freqbins-maxi*harms(i);
            freq_window=slice(minimum:maximum);
            [Pow,Ind]=max(freq_window);  
            % Interpolation--subtraction of 0.5 because we care about the
            % central frequency of the bin, not the upper boundary.
            if Ind==1 || Ind==length(freq_window) % edges
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
        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest); % weighted average of all harms.
    end
    pitch_data(:,found_note)=freqbin_estimate*(Nyquist/(freqbins-1));
end
