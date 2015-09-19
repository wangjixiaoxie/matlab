function pj_PSTH_Intan(dat,aux,datChan,auxChan,fs,stimLength,trialLength,binSizes,xNoise,saveFigs,plotRawDat,plotPSTH,plotRaster)
%% Make a PSTH and raster with user input data. For neural data recorded with the Intan system. 
% stimLength, trialLength, binSizes in ms!


dat = dat(datChan,:)';
aux = aux(auxChan,:)';


%% Get spike times (in seconds)
tL = trialLength;
sL = stimLength;
trialLength = trialLength*1e-3;
stimLength = stimLength*1e-3;
[pks,T,spikeThresh,timeVec] = pj_getSpikeTimes2(dat,fs,xNoise,plotRawDat,aux);

if saveFigs == 1
    h = gcf;
    t = ['','_rawData_detectedSpikes_channel',num2str(datChan),'_',num2str(spikeThresh),'spikeThreshold_withLaser.fig'];
    saveas(h, t)
end

%% Define windows around centered stimuli

C = midcross(aux,fs)
halfWay = (trialLength - stimLength)/2;
winBounds = C(1) - halfWay;

tSpikesPS = []; %Spike times with respect to stimulus window
for s = 1:(length(C)/2)
    winBounds = [winBounds, winBounds(end)+trialLength];
    newSpikes = T(T<winBounds(end) & T>winBounds(end-1));
    newSpikes = newSpikes-winBounds(end-1);
    tSpikesPS = [tSpikesPS;newSpikes];
end

numTrials = length(C)/2;

%% Plot PSTHs with given bin sizes (in ms)
if plotPSTH == 1
    for b = 1:length(binSizes)
        figure, hold on
        nBins = 1e3*trialLength/binSizes(b);
        intTest = ~mod(nBins,1);
        if intTest == 0
            disp('Warning: last bin is fractional!')
        end
        hist(tSpikesPS,nBins)
        X = [halfWay,(halfWay+stimLength)];
        Y = [0.01 0.01];
        line(X,Y,'Color','g','LineWidth',10)    % Draw line to indicate stimulus
        xlabel('Time (s)')
        ylabel('Spikes / Bin')
        title(sprintf(['PSTH for file ','',', spike threshold = ',num2str(spikeThresh),'\n'...
        'stim length = ',num2str(sL),'ms, trial period = ',num2str(tL),'ms, ',...
        'bin size = ',num2str(binSizes(b)),'ms\n',num2str(numTrials),' trials, laser head power = bla',num2str(datChan)]),...
        'FontWeight','bold')
        if saveFigs == 1
            h = gcf;
            t = [shortFn,'_PSTH_channel',num2str(datChan),'_',num2str(spikeThresh),'spikeThreshold_',num2str(binSizes(b)),'msBin.fig'];
            saveas(h, t)
        end
    end
end        

%% Plot raster of all detected spikes
if plotRaster == 1
    figure, hold on
    timePerLine = 10; %seconds
    sampsPerLine = timePerLine*fs;
    numLines = floor(length(timeVec)/sampsPerLine); %Very end of file not plotted
    nAux = aux/max(aux);
    for i = 1:numLines
        plot(timeVec(1:sampsPerLine),nAux(1:sampsPerLine),'g')
        oneLine = T(T<timePerLine);
        if i == 1
            yScat = ones(1,length(oneLine))*(i-1+1.5);
        else
            yScat = ones(1,length(oneLine))*(2*(i-1)+1.5);
        end
        scatter(T(T<timePerLine),yScat,2,'k.')
        nAux = nAux + 2;
        nAux = nAux(sampsPerLine:end);
        timeVec = timeVec - timePerLine;
        timeVec = timeVec(sampsPerLine:end);
        T = T(T>timePerLine);
        T = T - timePerLine;    
    end
    title(sprintf(['Raster plot for file yourfilehere, spike threshold = ',num2str(spikeThresh),'\n'...
        'stim length = ',num2str(sL),'ms, trial period = ',num2str(tL),'ms\n',...
        num2str(numTrials),' trials, laser head power = 30mW\n Right RA: 2.7mm Ventral, Channel ',num2str(datChan)]),...
        'FontWeight','bold')
    xlim([-5,timePerLine+5])
    ylim([-2,2*numLines+4])
    if saveFigs == 1
        h = gcf;
        t = [shortFn,'_Raster_channel',num2str(datChan),'_',num2str(spikeThresh),'spikeThreshold_',num2str(binSizes(b)),'msBin.fig'];
        saveas(h, t)
    end
end

%% Plot distribution of spike amplitudes
figure, hold on
%title(['Histogram of spike amplitudes'],'FontWeight','bold')
hist(pks,20)



end

