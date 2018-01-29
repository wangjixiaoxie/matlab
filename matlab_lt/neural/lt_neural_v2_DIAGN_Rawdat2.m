function lt_neural_v2_DIAGN_Rawdat(SummaryStruct, displaymode, skipnum)
%% lt 1/24/18 - 

%% lt 9/24/17 - for checking clustering. overlays spks with raw filtered dat.

% displaymode
%     order; % from song 1 to end
%     rand; % random, without replacement
%     skip; % from 1, skip N songs


%% MODIFY TO ITERATE
PlotSecondChan = 1;
plotcols={'m', 'r','c', 'b', 'g'};
i=1;
ii=1;
songnum = 5;


%% extract stuff [general]
cd(SummaryStruct.birds(i).neurons(ii).dirname);
birdname = SummaryStruct.birds(i).birdname;
exptname = SummaryStruct.birds(i).neurons(ii).exptID;

% -- collect spike times, song filenames, etc
load('times_data.mat'); % clusterclass
load('MetaDat.mat'); % metaDat

    cd ..
%% === list of songs, in order
songlist = 1:length(metaDat);
if strcmp(displaymode, 'rand')
    songlist = songlist(randperm(length(songlist)));
elseif strcmp(displaymode, 'skip')
    songlist = 1:skipnum:length(metaDat);
end

for songnum = songlist
    %%
    % -- load audio and raw neural for this song
    
    [amplifier_data,~,frequency_parameters, board_adc_data, ...
        ~ , amplifier_channels, ~] = pj_readIntanNoGui(metaDat(songnum).filename);
    indsamp=[amplifier_channels.chip_channel]==SummaryStruct.birds(i).neurons(ii).channel;
    
    neurdat = amplifier_data(indsamp, :);
    songdat = board_adc_data(1,:); clear board_adc_data;
    fs = metaDat(1).fs;
    X = (1:length(neurdat))/fs;
    
    %% plots
    hfig = lt_figure; hold on;
    hsplots = [];
    
    % ========== 1) plot audio
    hsplot = lt_subplot(6,1,1); hold on;
    title({[birdname '-' exptname '-n' num2str(ii)], [SummaryStruct.birds(i).neurons(ii).dirname(20:end)]}, 'Color','b');
    hsplots = [hsplots hsplot];
    plot(X, songdat);
    
    
    
    % ============= 2) plot neural
    hsplot = lt_subplot(6,1,2:3); hold on;
    title(['ch' num2str(SummaryStruct.birds(i).neurons(ii).channel) ...
        '-clust' num2str(SummaryStruct.birds(i).neurons(ii).clustnum)], 'Color','b');
    hsplots = [hsplots hsplot];
    [datfilt] =lt_neural_filter(neurdat, frequency_parameters.amplifier_sample_rate);
    plot(X, datfilt, 'k');
    
    % ============== overlay spike times
    % 1) boundaries of song, in cumulative ms
    cumsumsamps = [0 cumsum([metaDat.numSamps])];
    startms = (1000/fs)*cumsumsamps(songnum)+1;
    endms = (1000/fs)*cumsumsamps(songnum+1);
    
    % 2) extract spikes
    clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
    inds=find(cluster_class(:,1)==clustnum & ...
        cluster_class(:,2)>=startms ...
        & cluster_class(:,2)<=endms);
    
    for iii=1:length(inds)
        j=inds(iii);
        spktime=cluster_class(j, 2); % in ms
        spktime=spktime-startms+1;
        line([spktime spktime]/1000, [0 40]-i*20, 'Color', plotcols{i}, 'LineWidth',2);
    end
    
    % ================= NEURAL FOR OTHER CHANNELS
    hsplot = lt_subplot(6,1,4:5); hold on;
    hsplots = [hsplots hsplot];
    title([metaDat(songnum).filename ' (#' num2str(songnum) '/' num2str(length(metaDat)) ')'], 'Color','b');
    if PlotSecondChan==1
        count = 0;
        yrange = 500;
        numchans = length([amplifier_channels.chip_channel]);
        for j=1:numchans;
            
            chan2 = amplifier_channels(j).chip_channel;
            if chan2 == SummaryStruct.birds(i).neurons(ii).channel
                continue
            end
            
            count = count+1;
            [datfilt] =lt_neural_filter(amplifier_data(j, :), frequency_parameters.amplifier_sample_rate);
            plot(X, (datfilt-yrange)+yrange/(0.5*(numchans-1))+(count-1)*yrange/2, 'Color', [0.4 0.6 0.5]);
            lt_plot_text(mean(X), (-yrange)+yrange/(0.5*(numchans-1))+(count-1)*yrange/2, ['ch' num2str(chan2)], 'k')
            
        end
        axis tight;
    end
    
    linkaxes(hsplots, 'x');
    
    % ============= 3) timeline of all songs, showing where this song is
    lt_subplot(6,1,6); hold on;
    
    xx = [metaDat.song_datenum];
        xx = lt_convert_EventTimes_to_RelTimes(datestr(floor(xx), 'ddmmmyyyy'), xx);
        xx = xx.FinalValue;
    
    plot(xx, 1, 'ok');
    lt_plot(xx(songnum), 1, {'Color','r'});
    
    xlabel('time (day, rel first day)');
    
    axis tight
    ylim([0.5 1.5]);
    
        
    
    
    
    
    disp('paused ---');
    pause
    close(hfig);
end
