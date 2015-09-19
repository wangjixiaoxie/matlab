%% LT 7/31/14 - Modified:
% 1) plots a range of channgels (possibly) all, instead of one
%     (only supports one digital channel as of now)
    


function pj_PSTH_Intan(amplifier_data,board_dig_in_data,datChannels,digChan,fs,stimLength_ms,trialLength_ms,binSizes,xNoise,saveFigs,plotRawDat,plotPSTH,plotRaster)
% INPUTS:
% datChannels = [3, 5:10] (neural channels to look at);


%% Make a PSTH and raster with user input data. For neural data recorded with the Intan system. 
% stimLength, trialLength, binSizes in ms!

datChannels=sort(datChannels); 
num_DatChan=size(datChannels,2);
aux = board_dig_in_data(digChan,:)';
tL = trialLength_ms;
sL = stimLength_ms;
trialLength_s = trialLength_ms*1e-3;
stimLength_s = stimLength_ms*1e-3;

% Windows for PSTH
C = midcross(aux,fs); % tells where stimuli rise and fall
halfWay = (trialLength_s - stimLength_s)/2; % halfway between end of one and beginning of another stimulus.



for i=1:length(datChannels);
    datChan=datChannels(i);
    dat = amplifier_data(datChan,:)';



%% Get spike times (in seconds)
[pks{i},T{i},spikeThresh{i},timeVec{i}] = pj_getSpikeTimes2(dat,fs,xNoise,plotRawDat,aux); % timeVec - vector with values in seconds for each sample.


% if saveFigs == 1
%     h = gcf;
%     t = ['','_rawData_detectedSpikes_channel',num2str(datChan),'_',num2str(spikeThresh{i}),'spikeThreshold_withLaser.fig'];
%     saveas(h, t)
% end


%% Define windows around centered stimuli - CHANGE TO PLOT RASTERS, and not require inputs of durations. 


% get stim epoch boundaries (including baseline)
tSpikesPS = []; %Spike times with respect to stimulus window
winBounds = C(1) - halfWay; % Initiate winbounds
for s = 1:(length(C)/2) % for every stim epoch
    winBounds = [winBounds, winBounds(end)+trialLength_s];
    newSpikes = T{i}(T{i}<winBounds(end) & T{i}>winBounds(end-1));
    newSpikes = newSpikes-winBounds(end-1);
    tSpikesPS = [tSpikesPS;newSpikes];
end
tSpikesPS_out{i}=tSpikesPS; % CHANGE to separate all stim epochs, instead of concatenating.
numTrials = length(C)/2;

end


%% Plot PSTHs with given bin sizes (in ms)

PlotsPerFig=9;
[num_figs, numRows, numCols]=lt_get_subplot_size(num_DatChan,PlotsPerFig);

% Initiate figures;
for jj=1:num_figs;
    hfig_psth(jj)=figure; hold on; 
end

% Start plotting
c=1;
for i=1:length(datChannels);
    for b = 1:length(binSizes);
        figNum=ceil(c/PlotsPerFig);
        figure(hfig_psth(figNum));
        subplot(numRows(figNum), numCols(figNum),c-(figNum-1)*9);
        nBins = 1e3*trialLength_s/binSizes(b);
        % Make sure bin size conforms to trial lengths
        intTest = ~mod(nBins,1);
        if intTest == 0
            disp('Warning: last bin is fractional!')
        end
        % Plot histogram
        hist(tSpikesPS_out{i},nBins)
        X = [halfWay,(halfWay+stimLength_s)];
        Y = [0.01 0.01];
        line(X,Y,'Color','g','LineWidth',10)    % Draw line to indicate stimulus
        xlabel('Time (s)')
        ylabel('Spikes / Bin')
        %             title(sprintf(['PSTH for file ','',', spike threshold = ',num2str(spikeThresh{i}),'\n'...
        %                 'stim length = ',num2str(sL),'ms, trial period = ',num2str(tL),'ms, ',...
        %                 'bin size = ',num2str(binSizes(b)),'ms\n',num2str(numTrials),' trials, laser head power = bla',num2str(datChannels(i))]),...
        %                 'FontWeight','bold')
        title(['Channel ' num2str(datChannels(i)-1)]); % channel names go from 0 to 31
        c=c+1;
        %             if saveFigs == 1
        %                 h = gcf;
        %                 t = [shortFn,'_PSTH_channel',num2str(datChannels{i}),'_',num2str(spikeThresh{i}),'spikeThreshold_',num2str(binSizes(b)),'msBin.fig'];
        %                 saveas(h, t)
        %             end
    end
end
spikes_AllChans=[];
for i=1:size(tSpikesPS_out,2);
    spikes_AllChans=[spikes_AllChans; tSpikesPS_out{i}];
end

figure;
        hist(spikes_AllChans,nBins)
        X = [halfWay,(halfWay+stimLength_s)];
        Y = [0.01 0.01];
        line(X,Y,'Color','g','LineWidth',10)    % Draw line to indicate stimulus
        xlabel('Time (s)')
        ylabel('Spikes / Bin')
    title('all channels averaged');
    
    
%% Plot raster of all detected spikes
% if plotRaster == 1
%     figure, hold on
%     timePerLine = 10; %seconds
%     sampsPerLine = timePerLine*fs;
%     numLines = floor(length(timeVec)/sampsPerLine); %Very end of file not plotted
%     nAux = aux/max(aux);
%     for i = 1:numLines
%         plot(timeVec(1:sampsPerLine),nAux(1:sampsPerLine),'g')
%         oneLine = T(T<timePerLine);
%         if i == 1
%             yScat = ones(1,length(oneLine))*(i-1+1.5);
%         else
%             yScat = ones(1,length(oneLine))*(2*(i-1)+1.5);
%         end
%         scatter(T(T<timePerLine),yScat,2,'k.')
%         nAux = nAux + 2;
%         nAux = nAux(sampsPerLine:end);
%         timeVec = timeVec - timePerLine;
%         timeVec = timeVec(sampsPerLine:end);
%         T = T(T>timePerLine);
%         T = T - timePerLine;    
%     end
%     title(sprintf(['Raster plot for file yourfilehere, spike threshold = ',num2str(spikeThresh),'\n'...
%         'stim length = ',num2str(sL),'ms, trial period = ',num2str(tL),'ms\n',...
%         num2str(numTrials),' trials, laser head power = 30mW\n Right RA: 2.7mm Ventral, Channel ',num2str(datChan)]),...
%         'FontWeight','bold')
%     xlim([-5,timePerLine+5])
%     ylim([-2,2*numLines+4])
%     if saveFigs == 1
%         h = gcf;
%         t = [shortFn,'_Raster_channel',num2str(datChan),'_',num2str(spikeThresh),'spikeThreshold_',num2str(binSizes(b)),'msBin.fig'];
%         saveas(h, t)
%     end
% end

%% Plot distribution of spike amplitudes
% figure, hold on
% %title(['Histogram of spike amplitudes'],'FontWeight','bold')
% hist(pks,20)



end

