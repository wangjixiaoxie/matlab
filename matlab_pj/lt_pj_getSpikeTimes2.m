function [pks,T,spikeThresh,timeVec] = pj_getSpikeTimes2(dat,fs,xNoise,plotDetection,varargin)
%% Take raw neural recording, extract spike amplitudes and times.
%
% Inputs: dat = raw recording (vector); fs = sampling frequency (integer); xNoise = rms noise
% multiplier that determines spike threshold (float); 
% if plotDetection == 1, plots the raw data along with identified spikes;
% optionally also plots an auxiliary input, where the input is specified as
% the first and only argument of varargin (vector). 
%
% Outputs: pks = vector of identified spike amplitudes; T = corresponding vector of times at which pks occur (in ms);
% spikeThresh = threshold used to determine spikes; timeVec = vector of length(dat) with entry giving the sample number in ms. 

%% Set thresholds, calculate spike times.
noise = rms(dat(1:fs/2)); % First 500ms of data
spikeThresh = noise*xNoise; % May be modified

[pks,T] = findpeaks(dat,'minpeakheight',spikeThresh,'minpeakdistance',25);
Tms = T./fs;

timeVec = linspace(1/fs,length(dat)/fs,length(dat));
snapTimes = interp1(timeVec,timeVec,Tms,'nearest');
T = snapTimes;

%% Plot detected spikes with raw wave form
if plotDetection == 1
    figure, hold on
    plot(timeVec,dat,'b')
    scatter(T,pks,'r*')
    xlabel('Time (s)')
    if length(varargin) > 0
        aux = varargin{1};
        aux = (aux/max(aux))*mean(pks);
        plot(timeVec,aux,'g')
    end
end


    
    


end

