function lt_neural_concatExplore_v2(batchf, ChansToPlot, PlotRectDat, PlotFiltDat)

%% note: collects spikes individually for each song (i.e. within song threshold)

spkthresh_mult=3; % times median std.
Spk_PosOrNeg=-1; % 1 or -1
durPreSpike=0.7; % ms
durPostSpike=1; % ms

% NOTE: ONLY WORKS IF CHECKING A SINGLE CHANNEL!!!


%% lt 8/17/16 - exploratory only - concats neural (all chans) + song
% does not save anything
%
% ChansToPlot=neural chans (chip chan)

    windowsize=0.03; % from -2sd to +2sd [for smoothed rectified]

%% put names of all files in batch into a cell array

filenames=lt_batchsong_NamesToCell(batchf);



%%  concatenate channels as vectors - save metadata information along the way

fs_all=[];
ampDat_all=cell(length(ChansToPlot),1); % chan x data
songDat_all=[];
transsamps=[];


SpkWaveforms={}; % {songnum}(spknum, waveform);
MedianNoiseCollect=[];

% ----- COLLECT DATA
for i=1:length(filenames)
    
    fname=filenames{i};
    
    [amplifier_data,~ ,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(fname);
    
    % -- collect song
    songDat_all=[songDat_all board_adc_data(1,:)];
    transsamps=[transsamps length(board_adc_data(1,:))];
    
    % --- collect all amp dat
    for j=1:length(ChansToPlot)
        chan=ChansToPlot(j);
        ind=[amplifier_channels.chip_channel]==chan;
        dattmp=amplifier_data(ind, :); % dat for this chan
        
        ampDat_all{j}=[ampDat_all{j} dattmp]; % concat
        
        % --- collect all fs
        fs_all=[fs_all frequency_parameters.amplifier_sample_rate];
    end
    
    

    % ================= ADDED - COLLECT SPIKES FOR THIS SONG
    % ------ what is num inds tha make up spike duration? (for extracting
% spike)
    SpkPreSamps=(durPreSpike/1000)*frequency_parameters.amplifier_sample_rate;
    SpkPostSamps=(durPostSpike/1000)*frequency_parameters.amplifier_sample_rate;

    % -- filter
    datfilt=lt_neural_filter(dattmp, frequency_parameters);
    
    % === get RMS
    MedianNoise=median(abs(datfilt)./0.6745); % median noise (not std), from http://www.scholarpedia.org/article/Spike_sorting
    SpikeThreshold=spkthresh_mult*MedianNoise;
    MedianNoiseCollect=[MedianNoiseCollect MedianNoise];
    
    [SpikePks, SpikeInds]=findpeaks(Spk_PosOrNeg*datfilt,'minpeakheight',...
        SpikeThreshold,'minpeakdistance',floor(0.0003*frequency_parameters.board_adc_sample_rate)); % 0.3ms min separation btw peaks.
    
    numsampstmp=length(datfilt);
    % --- collect spike waveforms
    count=1;
    for j=1:length(SpikeInds);
        
        if SpikeInds(j)>SpkPreSamps & SpikeInds(j)+SpkPostSamps-1<=numsampstmp
        spkdat=datfilt(SpikeInds(j)-SpkPreSamps:SpikeInds(j)+SpkPostSamps-1);
        SpkWaveforms{i}(count, :) = spkdat; % {songnum}(spknum, waveform)
        count=count+1;
        end
    end
end


% --- confirm that all have same sample rate
if length(fs_all)>1
    assert(all(diff(fs_all)==0), 'problem, FS not same for all file')
end

for i=1:length(filenames)
       disp(['Median noise for song num ' num2str(i) ' = ' num2str(MedianNoiseCollect(i))]);
 
end

disp(['== Sample rate = ' num2str(fs_all(1))]);


%% PLOT
figcount=1;
subplotrows=2*length(ChansToPlot)+1;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];

% - song, raw
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots=[hsplots hsplot];
tt=[1:length(songDat_all)]/fs_all(1);
plot(tt, songDat_all);

%  line for each transition
transsamps=cumsum(transsamps);
for i=1:length(transsamps)
    x=tt(transsamps(i));
    line([x x], ylim, 'Color', 'r')
end

% - plot each chan
for i=1:length(ChansToPlot)
    chan=ChansToPlot(i);
    
    % === 1) neural, filtered
    if PlotFiltDat==1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    title(['chan ' num2str(chan)]);
    
    datfilt=lt_neural_filter(ampDat_all{i}, frequency_parameters);
    plot(tt, datfilt, 'k');
    end
    
    % === 2) neural, smoothed
    if PlotRectDat==1
    % === rectify
    datfilt=abs(datfilt);
    
    % == smoooth
    % Construct a gaussian window
    sigma=(windowsize/4)*fs_all(1); %
    numsamps=4*sigma; % (get 2 std on each side)
    alpha= numsamps/(2*sigma); % N/2sigma
    gaussFilter = gausswin(numsamps, alpha);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    dat_smrect = conv(datfilt, gaussFilter);
    
    % == PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    plot(tt, dat_smrect(numsamps/2:end-numsamps/2), 'k');
    hsplots=[hsplots hsplot];
    end
end




linkaxes(hsplots, 'x');

%% --- PLOT SPIKE WAVEFORMS


% ==== V1, lines 
figcount=1;
subplotrows=3;
subplotcols=8;
fignums_alreadyused=[];
hfigs=[];


numrandspikes=100;
hsplots=[];

for i=1:length(SpkWaveforms)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    % ==== overlay subset of spikes
    inds=randperm(size(SpkWaveforms{i},1), numrandspikes);
    
    dattmp=SpkWaveforms{i}(inds, :);   
    plot(dattmp', 'b-');
%     plot(dattmp', 'Color', [0 0 0 0.5]);
    
    wave_mean=mean(dattmp);
    plot(1:length(wave_mean), wave_mean, 'k','LineWidth', 3);
    axis tight
end

linkaxes(hsplots, 'xy');
% subplotsqueeze(gcf, 1.3);


% ==== V2 -- dots, to see density
figcount=1;
subplotrows=3;
subplotcols=8;
fignums_alreadyused=[];
hfigs=[];


numrandspikes=300;
hsplots=[];

for i=1:length(SpkWaveforms)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    % ==== overlay subset of spikes
    inds=randperm(size(SpkWaveforms{i},1), numrandspikes);
    
    dattmp=SpkWaveforms{i}(inds, :);   
    plot(dattmp', 'b:');
%     plot(dattmp', 'Color', [0 0 0 0.5]);
    
    wave_mean=mean(dattmp);
    plot(1:length(wave_mean), wave_mean, 'k','LineWidth', 3);
    axis tight
end

linkaxes(hsplots, 'xy');
% subplotsqueeze(gcf, 1.3);


