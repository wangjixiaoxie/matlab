function lt_neural_CheckClust(varargin)
% e.g. lt_neural_CheckClust(11, 'Batch1438to1909', 14)
% RUn in day folder. must have lower directory properly named

% varargin{1} = channel
% varargin{2} = batchfile
% varargin{3} = 2nd channel, to check common noise artifact. [optional]

% HAND REMOVE SPIKES - CHANGES THEIR CLUSTER NUMBER TO -1, automatically
% saves one each movement to new song file


%% 
% NOTE: in wave-clus, any removed spikes will be put into cluster 0 (i.e
% noise cluster) (if you use "save clusters" in wave-clus) opening
% wave-clus does nothing. 
% Have verified that -1 is converted to 0.

% Have verified that indeed spikes are removed from the correct
% cluster.

% Have verified that spieks are removed and saved correctly. 

% Have verified that whether or not run (and save) wave_clus post
% hand-checking still works fine for analysis (will just change -1 to 0,
% clusnums)


%%

global mm
global spkinfo
global spikes_cat
global cumulativetime_ms

cumulativetime_ms = 0;
mm =1;

%% plot
channel_board = varargin{1};
batchf = varargin{2};

% --- go to data folder
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)

% ---- load, go back up, and plot.
% 1) EXTRACT AND CONCAT ALL SONG FILES and neural files
load('MetaDat.mat'); % contains file names
% --- load concat neural and spikes
neural_cat=load('data.mat');
spikes_cat=load('times_data.mat');

% --- SAVE A BACKUP VERSION IF NOT ALREADY
if ~exist('times_data_BackupBeforeHandCheck.mat', 'file')
eval('!cp times_data.mat times_data_BackupBeforeHandCheck.mat');
end

cd ..

lt_figure; hold on;
plot_song_file(varargin{1}, varargin{2}, varargin{3}, '', metaDat, neural_cat);
%% === put a ui button

hremovespikes = uicontrol('Style', 'pushbutton', ...
    'String', 'Remove Spikes', 'Position', [1 20 100 30], ...
    'Callback', @removespikes);

SecondChan = varargin{3};
Pos = get(gcf, 'Position');

axespos = get(gca, 'Position');
hnextsong = uicontrol('Style', 'pushbutton', ...
    'String', 'Next song (and save)', 'Position', [300 1 100 30], ...
    'Callback', {@nextsongcallback, channel_board, batchf, ...
    metaDat, neural_cat, SecondChan});


% %% ======
% while (1)
%     skipsong =0;
% while skipsong==0
%     
%     set(gcf, 'Units', 'Normalized')
%     disp('WAITING!!!')
%     waitforbuttonpress
%     finalrect = rbbox;
%     
% %     annotation('rectangle', finalrect);
%     
%     %===== convert from figure units to units in the figure
%     set(gca, 'Units', 'normalized')
%     axes_pos = get(gca, 'Position'); % position of axes rel to figure (norm)
%     
%     Xlim = get(gca, 'XLim'); % axes limits (actual units)
%     Ylim = get(gca, 'Ylim');
%     xrange = Xlim(2)-Xlim(1);
%     
%     % -- actual unit for left limit
%     tmp = (finalrect(1) - axes_pos(1))/axes_pos(3);
%     X1 = Xlim(1)+tmp*xrange
%     
%     plot(X1, 0.3, 'og');
%     
%     % -- actual unit for right limit
%     tmp = (finalrect(1)+finalrect(3) - axes_pos(1))/axes_pos(3);
%     X2 = Xlim(1) + tmp*xrange
%     
%     plot(X2, 0.3, 'or');
%     
% end
% 
% % === move to next file.
% close all; 
% lt_figure; hold on;
% % === put a ui button
% 
% axespos = get(gca, 'Position');
% hnextsong = uicontrol('Style', 'pushbutton', ...
%     'String', 'Next song', 'Position', [1 1 100 30], ...
%     'Callback', @nextsongcallback);
% 
% mm = mm +1;
% plot_song_file(varargin{1}, varargin{2}, '', mm, metaDat, neural_cat, ...
%         spikes_cat);
% 
% 
% % === turn off skip song
% skipsong = 0;
% end
end

function removespikes(source, eventdata)
global spkinfo 
global spikes_cat

    set(gcf, 'Units', 'Normalized')
    disp('WAITING!!!')
    waitforbuttonpress
    finalrect = rbbox;
    
%     annotation('rectangle', finalrect);
    
    %===== convert from figure units to units in the figure
    set(gca, 'Units', 'normalized')
    axes_pos = get(gca, 'Position'); % position of axes rel to figure (norm)
    
    Xlim = get(gca, 'XLim'); % axes limits (actual units)
    Ylim = get(gca, 'Ylim');
    xrange = Xlim(2)-Xlim(1);
    
    % -- actual unit for left limit
    tmp = (finalrect(1) - axes_pos(1))/axes_pos(3);
    X1 = Xlim(1)+tmp*xrange
    
    plot(X1, 0.3, 'og');
    
    % -- actual unit for right limit
    tmp = (finalrect(1)+finalrect(3) - axes_pos(1))/axes_pos(3);
    X2 = Xlim(1) + tmp*xrange
    
    plot(X2, 0.3, 'or');
    
    
    % ============ REMOVE ALL SPIKES WITH TIMEPOINT BETWEEN X1 AND X2
    % find all spikes (from all clusters) that fall within this window
    indstmp = spkinfo.spktimesRelFile_all > X1*1000 ...
        & spkinfo.spktimesRelFile_all < X2*1000;
    
    clustremove = spkinfo.cluster_all(indstmp);
    spktimesglobremove = (1/1000)*spkinfo.spktimesGlobal_all(indstmp);
    spktimesfileremove = (1/1000)*spkinfo.spktimesRelFile_all(indstmp);
    
    disp(['===== CHANGED SPIKES CLUST TO -1']);
    disp(['Times (within file):' num2str(spktimesfileremove')]);
    disp(['Times (global): ' num2str(spktimesglobremove')]); 
    disp(['Clusters : ' num2str(clustremove')]);
    
    
    % ----- ACTUALLY REMOVE, BY CHANGING CLUSTER TO -1
    indsGlobalSpksToRemove = spkinfo.inds_all(indstmp);
    spikes_cat.cluster_class(indsGlobalSpksToRemove,1) = -1;
    
end


function nextsongcallback(source, eventdata,  channel_board, batchf, metaDat, neural_cat, SecondChan)
global spikes_cat
global mm

% === claer
htmp = get(gcf, 'Children');
for i=1:length(htmp)
    cla(htmp(i));
end


% === save spikes_cat
% --- go to data folder
ver = '';
currentver = version;
if currentver(1) >= 7
    ver = '-v6';
end
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)
save('times_data.mat', '-struct', 'spikes_cat' , ver);
disp('SAVED');
cd ..


% === load next song.
mm = mm +1;
plot_song_file(channel_board, batchf, SecondChan, '', metaDat, neural_cat, ...
        spikes_cat);

end



function plot_song_file(varargin)
global spkinfo
global spikes_cat
global mm
global cumulativetime_ms
% clear spkinfo % so can redetermine spkinfo.

set(gcf, 'toolbar', 'figure');

channel_board = varargin{1};
batchf = varargin{2};
% mm = varargin{4};
metaDat = varargin{5};
neural_cat = varargin{6};
% spikes_cat = varargin{7};

plotcols={'m', 'r','c', 'b', 'g'};

PlotSecondChan=0;
if ~isempty(varargin{3})
    PlotSecondChan=1;
    SecondChan = varargin{3};
end

% % --- go to data folder
% datdir=['Chan' num2str(channel_board) 'amp-' batchf];
% cd(datdir)
% 
% % ---- load, go back up, and plot.
% % 1) EXTRACT AND CONCAT ALL SONG FILES and neural files
% load('MetaDat.mat'); % contains file names
% % --- load concat neural and spikes
% neural_cat=load('data.mat');
% spikes_cat=load('times_data.mat');
% cd ..


AllSongs_old=[];
AllNeural_old=[];
AllLabels=[];
AllOnsets=[];
AllOffsets=[];
cumulative_filedur=0; % keeps track as concatenating

AllNeural_SecondChan = [];


% --- load audio and old neural (from actual files)
hsplots = [];
% -- load original sound and neural
[amplifier_data,~,frequency_parameters, board_adc_data, ...
    board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(mm).filename);

indsamp=find([amplifier_channels.chip_channel]==channel_board);
neurdat = amplifier_data(indsamp, :);

% -- load labels, onsets, offsets
tmp=load([metaDat(mm).filename '.not.mat']);
AllLabels=[AllLabels tmp.labels];

tt = [1:length(board_adc_data(1,:))]/frequency_parameters.amplifier_sample_rate;

datdur_ms = 1000*length(board_adc_data(1,:))/frequency_parameters.amplifier_sample_rate;

% - a. song raw
hsplot = lt_subplot(3,1,1); hold on;
plot(tt, board_adc_data(1,:));
hsplots=[hsplots hsplot];

% onsets and labels
for i=1:length(tmp.onsets)
    line([tmp.onsets(i) tmp.offsets(i)]/1000, [1.65 1.65], 'LineWidth', 2, 'Color', 'r');
    lt_plot_text(tmp.onsets(i)/1000, 1.5, tmp.labels(i), 'r')
end


if PlotSecondChan==1
    ind =find([amplifier_channels.chip_channel]==SecondChan);
    hsplot = lt_subplot(3,1,2); hold on;
    [datfilt] =lt_neural_filter(amplifier_data(ind, :), frequency_parameters.amplifier_sample_rate);
    plot(tt, datfilt, 'k');
    hsplots=[hsplots hsplot];
    title(['second chan (' num2str(SecondChan) ')']);
    
end


% - d. neural (cat)
% hsplot = lt_subplot(2,1,2); hold on;
hsplot = lt_subplot(3,1,3); hold on;
[datfilt] =lt_neural_filter(neurdat, frequency_parameters.amplifier_sample_rate);
plot(tt, datfilt, 'k');
hsplots=[hsplots hsplot];


% overlay spike times (cat)
numclust=max(spikes_cat.cluster_class(:,1));

inds_all = []; % ----- SAVE INFORMATION, IN CASE WANT TO THROW OUT SPIKES
spktimesRelFile_all = [];
spktimesGlobal_all = [];
cluster_all = [];
for i=1:numclust
    inds=find([spikes_cat.cluster_class(:,1)]==i & ...
        spikes_cat.cluster_class(:,2)>cumulativetime_ms ...
        & spikes_cat.cluster_class(:,2)<=cumulativetime_ms+datdur_ms);
    
    inds_all = [inds_all; inds];
    cluster_all = [cluster_all; ones(length(inds),1)*i];
    for iii=1:length(inds)
        ii=inds(iii);
        spktime=spikes_cat.cluster_class(ii, 2); % in ms
        spktimesGlobal_all = [spktimesGlobal_all; spktime];
        
        spktime=spktime-cumulativetime_ms;
        line([spktime spktime]/1000, [0 40]-i*20, 'Color', plotcols{i}, 'LineWidth',2);
        spktimesRelFile_all = [spktimesRelFile_all; spktime];
    end
end
spkinfo.inds_all = inds_all;
spkinfo.spktimesGlobal_all = spktimesGlobal_all;
spkinfo.spktimesRelFile_all = spktimesRelFile_all;
spkinfo.cluster_all = cluster_all;


% --- link all
linkaxes(hsplots, 'x');

%         if PlotSecondChan==1
%             ind=find([amplifier_channels.chip_channel]==SecondChan);
%             AllNeural_SecondChan=[AllNeural_SecondChan amplifier_data(ind, :)];
%         end
cumulativetime_ms = cumulativetime_ms + datdur_ms;

lt_subtitle(['song #' num2str(mm) '/' num2str(length(metaDat)) '(' metaDat(mm).filename ')'])
% --------------------------
end

%% TROUBLESHOOTING
% 1) IF want to undo, need to go out of gui, load time_data into tmp
if (0)
Numsongs = 62;
cumtime = 0;
for i=1:Numsongs
    
    cumtime = cumtime + metaDat(i).numSamps/30000;
    
end


inds = find(tmp.cluster_class(:,2)/1000 > 2497.3 ...
    & tmp.cluster_class(:,2)/1000 < 2529.312 ...
    & tmp.cluster_class(:,1) == -1);

tmp.cluster_class(212148:212157, 1)=1;
end
