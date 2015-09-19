function pj_antidromicIntan(file,stimChan,winSize,varargin)
%% Take an intan .rhd file as input, detect stim (digital input on stimChan),
% average recording on all channels around the artifact, plot in a map
% corresponding to the location of the sites on the probe. Useful for
% antidromic stimulation. 

%% Specify mapping from recording channel to probe site
vertMap = [32,4,28,20,18,14,6,2, ...
           24,22,10,12,16,26,8,30, ...
           31,13,27,21,7,15,9,1, ...
           25,5,17,11,3,25,19,29];

horzMap = [16,2,14,10,9,7,3,1, ...
           12,11,5,6,8,13,4,15, ...
           32,23,30,27,20,24,21,17, ...
           28,19,25,22,18,29,26,31];


[ampDat,digDat,freqParams] = pj_readIntanNoGui(file);
fs = freqParams.amplifier_sample_rate;

if length(varargin) > 0
    if varargin{1} == 'bandpass'
        for c = 1:32
            ampDat(c,:) = bandpass_filtfilt(ampDat(c,:),fs,300,10000);
        end
    end
end

stim = digDat(stimChan,:);
cross = midcross(stim);
cross = cross(1:2:end);
cross = floor(cross);

winSamps = fs*1e-3*winSize;
windowed = NaN(32,winSamps,length(cross)); %winSize = ms around each stim

for d = 1:size(ampDat,1)
    for c = 1:length(cross)
        winStart = cross(c) - winSamps/2;
        windowed(d,:,c) = ampDat(d,winStart:winStart+winSamps-1);
    end
end

meanWin = mean(windowed,3); %Average across stimulations
figure, plot(meanWin(1,:)) %Does result for one channel look reasonable?
figure, plot(ampDat(3,:))
figure, plot(ampDat(5,:))

timeVec = linspace(1e3/fs,1e3*(winSamps/fs),winSamps); %ms

%figure
%ax = [];
%for d = 1:size(ampDat,1)
%    ax(d) =  subplot(2,16,horzMap(d));
%    plot(timeVec,meanWin(d,:));
%    xlabel('Time (ms)')
%    ylabel('Voltage (uV)')
%    linkaxes(ax);
%end

for j = 1:5
    figure, hold on
    for i = 1:size(windowed,3)
        plot(timeVec,windowed(j,:,i))
        xlabel('Time (ms)')
        ylabel('Voltage (uV)')
        ylim([-500 500])
    end
end
    
    
    







end

