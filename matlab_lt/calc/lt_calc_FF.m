function [FF, PC, T] = lt_calc_FF(syldat, fs, Fwind, Twind, plotSpec)

% INPUTS -----
% syldat % acoustic syl dat (raw, unfiltered is fine)
% fs % sample rate
% Fwind % [F_low F_high] frequency window within which to find FF
% Twind % time window to take average PC [mintime maxtime], seconds

% OUTPUTS---
% FF % pitch
% PC % vector, pitch controu
% t % vector, time base for PC
%

if ~exist('plotSpec', 'var')
    plotSpec = 0;
end
    

N=1024;
sigma=1;
SAMPLING=fs;
OVERLAP=1020;

F_low=Fwind(1);
F_high=Fwind(2);

mintime=Twind(1);
maxtime=Twind(2);

%%


t=-N/2+1:N/2;
sigma=(sigma/1000)*SAMPLING;
%Gaussian and first derivative as windows.
w=exp(-(t/sigma).^2);
timebins=floor((length(syldat)-(N))/(N-OVERLAP))+1;  % number of bins, i.e. how many times the window N will be shifted to reach the end of the data.
freqbins=(N/2)+1;

Nyquist=SAMPLING/2;
step=Nyquist/freqbins; % This is the width of freqency bins (in Hz)
mini=round(F_high/step);
maxi=round(F_low/step); % step numbers for the min and max of the window you are interested in.

sonogram=zeros(freqbins,timebins);
[S,F,T]=spectrogram(syldat, w,OVERLAP,N,SAMPLING);
sonogram=abs(flip(S,1)); % gaussian windowed spectrogram, for entire data duration
F=flip(F);

% -
freqbin_estimate=nan(1, size(sonogram,2));
for currenttime_bin=1:size(sonogram,2)
    slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
    %         Powerful(1)=0;
    minimum=freqbins-mini;
    maximum=freqbins-maxi;
    freq_window=slice(minimum:maximum);
    
    [Pow,Ind]=max(freq_window);
    
    % Interpolation--subtraction of 0.5 because we care about the
    % central frequency of the bin, not the upper boundary.
    if Ind==1 || Ind==length(freq_window) % edges
        Indest=Ind;
    else
        Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
    end
    
    %                 disp(['central freq bin: ' num2str(Ind) ' - ' num2str(Indest)]);
    
    Real_Index=minimum+Indest-1; % -1 to account for window;
    Freqbinest=freqbins-Real_Index;
    %             Powerful=Pow;
    %         normalizer=sum(Powerful);
    %         normpower=Powerful/normalizer;
    %         freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest); % weighted average of all harms.
    
    freqbin_estimate(currenttime_bin)=Freqbinest; % weighted average of all harms.
end
PC=freqbin_estimate*(Nyquist/(freqbins-1));
PC=single(PC); % less memory use

% ==== extract FF
[~, minind]=min(abs(T-mintime));
[~, maxind]=min(abs(T-maxtime));

FF=mean(PC(minind:maxind));
% ----------------------------

% === DEBUG - plot spectrogram, with PC overlayed
if plotSpec==1
    %                 if strcmp(syl, 'b') & length(PC)==219
    %                 if counter==1781;
    lt_figure; hold on;
    %                     title(syl);
    imagesc(T, F, log(sonogram)); hold on;
    plot(T, PC, 'k', 'LineWidth', 2);
    line([mintime mintime], ylim, 'Color', 'm');
    line([maxtime maxtime], ylim, 'Color', 'm');
    line([mintime maxtime], [FF, FF], 'Color','r');
    
    line([mintime maxtime], [F_low F_low], 'Color', 'm');
    line([mintime maxtime], [F_high F_high], 'Color', 'm');
    axis tight;
    
    pause; close all;
    %                 end
end
% OUTPUT: 4 things