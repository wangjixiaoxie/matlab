%% TO DO:
% 1) save new time base for amplitude plotting (lin warped)
% 2) 


%% ===== plot song

songname='br92br54_170416_111151.rhd';
% songnotmat='bk7_160809_110736.rhd.not.mat';

lt_neural_extractaudio(songname);

% ---- overlay onsets and offsets [first label with evsonganaly]
notmatstruct=load(songnotmat);

for i=1:length(notmatstruct.onsets)
    
    line([notmatstruct.onsets(i) notmatstruct.offsets(i)]./1000, [1.65 1.65], 'Color','r','LineWidth', 7);
end


%% get a list of sampling rates of all files


%% ====== plot single file dat [align neural and song]
close all;
filename='pu69wh78_171120_205155.rhd';
ChansToPlot.DigChans_zero=[0]; % make string "all" to plot all that exist. empty array to ignore
ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
% ChansToPlot.AmpChans_zero=[9 14 19];
ChansToPlot.AmpChans_zero=[9 11 12 14 1 18];
% ChansToPlot.AmpChans_zero=[8 9 11 16 17 20 21];
ChansToPlot.AmpChans_zero=[9 14 17 18 21];
ChansToPlot.AmpChans_zero=[9 14 21];
ChansToPlot.AmpChans_zero=0:31;

% neuralFiltLow=500;
neuralFiltLow=300;

PlotWhat.raw=0;
PlotWhat.filt=1;
PlotWhat.rect_sm=0;
PlotWhat.raster=0;
PlotWhat.digital=0;

Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
Raster.ThrXNoise=4; % threshold for units, used for all channels, uses absolute data for peak detection
Raster.PosOrNeg=-1; % -1 or +1, for crossings.

        numsubplots = 3;

lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster, ...
            numsubplots)

%% ======= clean song files




%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ANALYSIS PIPELINE
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ==================================================== VARIOUS THINGS
%% PCA to try to remove common noise

% - SAVES ONE CHANNEL, READY FOR WAVE_CLUS (HERE SAVES ALREADY BANDPASS
% FILTERED, WHILE WAVECLUS ITSELF ALSO RUNS A FILTER)
% - DENOISES THAT CHANNEL CASED ON LINEAR REGRESSION AGAINST THE OTHER
% CHANNELS + CHANNEL MEDIAN (TO IGNORE SPIKES - I.E. IF NOISE THEN ALL
% CHANNELS SHOULD BE DOING SAME THING)


close all;
% ChansToUse=[9 15 18 11 14];
ChansToUse=[9 15 18 11];
batchf='Batch1710to1925';
plotOn = 1;
PlotAllChanSep = 0; % then makes one plot for each chan, comparing diff methods

onlyPlotFinalRegressionFig = 1; % overwrites other plot params, only plots final summary fig.
skipPCAmethod = 1; % then only does regression method - useful if just want to run fast and save
save_raw_dat = 1;
ChanToSave = 9;
Subsample_NumSec = 160; % number of seconds to use for building model.

lt_neural_commonNoise(batchf, ChansToUse, plotOn, save_raw_dat, ...
    PlotAllChanSep, ChanToSave, onlyPlotFinalRegressionFig, skipPCAmethod, Subsample_NumSec)

% NOTE (12/2/16) - currently regression method worked well for one set of
% songs and PCA method for another. Neither method worked great. 

% NOTE (2/22/17) - regression method seems to work well. An issue that I am
% worried about is introducing spikes where there should not be any. Is not
% easy to confirm that this method is not producing artifacts. So for now
% will avoid using this. But by eye looks like dramatically reduces one of
% the sources of noise, but there is another source of noise in many files
% that is not removed, and is actually enhanced in some cases - i.e.
% different weights depending on source.



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==================================================== PREPROCESSING
clear all; close all;
% channel_board = [8 11 14 18 20]; % wh6
% channel_board = [9 11 12 14 15 18]; % bu77
% channel_board = [14];
channel_board = [9 14 17 18 21];
channel_board = 0:31;
channel_board = 18;
batchf = 'Batch1244to2119';

%% ==== exploratory - concat all audio and neural and plot for each neural channel
close all;

% -----  v1, plots just raw neural, filtered
% batchtoplot=batchf;
if (0)
    lt_neural_concatExplore(batchf, channel_board);
end

% ----- v2, plots filtered neural, smoothed, and spike waveforms
PlotRectDat=0; % 1, plots, 0 skips.
PlotFiltDat=1; % usually 1, filt neural.
PosAndNeg =0; % then gets both. if 0, then just downwards


lt_neural_concatExplore_v2(batchf, channel_board, PlotRectDat, PlotFiltDat, PosAndNeg); % WAVEFORMS ONLY PLOTTED FOR ONE CHANNEL!!

% ------- v3 -plots song by song
if (0)
    fid = fopen(batchf);
    filename = fgetl(fid);
    while ischar(filename)
        
        
        ChansToPlot.DigChans_zero=[0]; % make string "all" to plot all that exist. empty array to ignore
        ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
        % ChansToPlot.AmpChans_zero=[9 14 19];
        ChansToPlot.AmpChans_zero=channel_board;
        % ChansToPlot.AmpChans_zero=[8 9 11 16 17 20 21];
        % ChansToPlot.AmpChans_zero=[8 9 16 17 20 21];
        
        % neuralFiltLow=500;
        neuralFiltLow=300;
        
        PlotWhat.raw=0;
        PlotWhat.filt=1;
        PlotWhat.rect_sm=0;
        PlotWhat.raster=0;
        PlotWhat.digital=0;
        
        Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
        Raster.ThrXNoise=3; % threshold for units, used for all channels, uses absolute data for peak detection
        Raster.PosOrNeg=-1; % -1 or +1, for crossings.
        
        numsubplots = 4;
        
        lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster, numsubplots)
        
        disp([' ======== ' filename]);
        
        pause;
        filename = fgetl(fid);
    end
end


%% ==== concatenate multiple files for given channel, [and saves]
% based on expectation of duration for single unit.
% -- saves into downstream folder
close all;
lt_neural_concatOneChan(batchf, channel_board)
cd(['Chan' num2str(channel_board) 'amp-' batchf]);


%% ==== run wave_clus on this concatted data



%% ======= CREATE NOT.MAT. FILE FOR ALL SONG FILES IN BATCH [AUTO]

lt_neural_AutoMakeNotmat(batchf);

        
%% ==== [SANITY CHECK] plot song files in entirety, aligned to extracted spikes
% -- makes multiple plots if too much dat.
close all;
PlotSecondChan = 1;
SecondChan = 17;
plotcols={'m', 'r','c', 'b', 'g'};

% want to plot 2nd channel to compare noise?
if (0)
    lt_neural_AlgnWavclus(batchf, channel_board, plotcols, PlotSecondChan, SecondChan);
end
        
% === VERSION 2 - PLOTS EACH SONG FILE ONE AT A TIME - CAN CLOSE THEM ONE
% AT A TIME
maxfiguuresopen = 20;
figsstepsize = 5; % num figs to jump thru, useful if many song. 1 is default. DOES NOT WORK
lt_neural_AlgnWavclus_v2(batchf, channel_board, plotcols, PlotSecondChan, SecondChan, maxfiguuresopen, figsstepsize);


%% ===== HAND REMOVE NOISY SPIKES
% IMPORTANT!!: only run this once you are satisfied with wave-clus cluster
% assignments. If reassign clusters in wave_clus, then will delete all
% hand-removed clusters (cluster -1)
% NOTE: can check which spikes removed using lt_neural_AlgnWavclus_v2

close all;

SecondChan = 18;
lt_neural_CheckClust(channel_board, batchf, SecondChan);


% ========== TROUBLESHOOTING
% i.e. made mistake of hand checking before finalized (consolidated)
% clusters. THen run below (MODIFIED!) to manually change cluster numbers,
% whithout wiping out hand checked (cluster -1)

if (0) 
tmp = load('times_data_TMPBACKUP.mat');

% to convert cluster 2 to 0
inds = tmp.cluster_class(:,1)==2;
tmp.cluster_class(inds,1) = 0;

% to convert clusters 3 and 4 to 1
inds = tmp.cluster_class(:,1) == 3 | tmp.cluster_class(:,1) == 4;
tmp.cluster_class(inds,1) = 1;

% to save
ver = '-v6';
save('times_data.mat', '-struct', 'tmp' , ver);
disp('SAVED');
end
%% +++++++++++++++++++++++++++++++++++++== NOTES THINGS TO UPDATE WITH ANALYSIS PIPELINE

% save all things to a common folder. Levels:
% bird, 
% batch - metadata (real times of files. names of files, etc).
%     do not save neural data. [if eventually need this can easily write
%     script to get]
%     do not save song data
% labels. onsets, offsets. [write script to update all labels]. or just use
% labels from notmat files
% spike times
% PROCESSED DATA: e.g. pitch ...




%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ++++++++++ ANALYSIS (NEURONS ONE BY ONE) 
%% ==== EXTRACT SONG, LABEL,SongDat ONSETS, SPIKE DATA
close all; clear all;
batchf='Batch0922to1251';
channel_board=14;
extractsound = 1;
[SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, extractsound);


%% ===== Extract song segments that match the desired regexp 
% along with associated neural data
close all;

% - desired motifs
% regexpr_str='[^vb](b)'; % token determines where align. predur is relative to token.
% regexpr_str='(g)h'; % token determines where align. predur is relative to token.
% regexpr_str='n(h)h'; % token determines where align. predur is relative to token.
% regexpr_str='[^h](h)h'; % token determines where align. predur is relative to token.
regexpr_str='nk(h)'; % token determines where align. predur is relative to token.
predur=0.25;
postdur=0.1;
alignByOnset=1;

% - entire motifs
% regexpr_str='(S)[\w-]*?E'; % gets any S ---- E
% predur=6; % sec
% postdur=8; % sec
% alignByOnset=1;

% - Extract song motifs (onset and offset) automatically)
% regexpr_str='WHOLEBOUTS';
% predur=8; % sec
% postdur=3; % sec
% alignByOnset=0;

keepRawSongDat = 1;
[SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
    regexpr_str, predur, postdur, alignByOnset, [], [], keepRawSongDat);



%% ===== [OPTIONAL] extract segments automatically - motif on and off


%% ===== uniformly stretch durations of segments
% assumes that onset of syl 1 to offset of syl last should be aligned.
close all;
TargetDur=[]; 
[SegmentsExtract, Params] = lt_neural_LinTimeWarp(SegmentsExtract, Params, TargetDur);


%% ==== plot rasters

close all;
useRescaled=0; % 1, then need to run LinTimeWarp first (plots scaled spikes, not song dat)
plotAllSegs=1; % then plots each trial own plot.
[Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs);


%% === compare rasters and mean firing rate for differnet motifs
% will align them by their own respective alignment locations.
% - can either rescale them to all be the same dur, or not rescale.

close all; 
% MotifList_regexp={'(S)[\w-]*?E'};
% MotifList_regexp={'g(h)h', 'n(h)h'};
% MotifList_regexp={'v(b)b', '[^vb](b)b'};
% MotifList_regexp={'n(b)b', 'v(b)b'};
% MotifList_regexp={'[nv](b)', '[gm](b)'};
% MotifList_regexp={'(v)', '(b)'};
% % MotifList_regexp={'[^b](b)', 'b(b)'};
% MotifList_regexp={'(g)b', '(g)h'};
% MotifList_regexp={'gh(h)hh'};
MotifList_regexp={'g(h)', 'n(h)', 'h(h)'};

predur=0.1;
postdur=0.1;
LinearWarp=0; % 0=no; 1=yes, each motif individually; 2=yes, all motifs to global motif median dur
suppressplots=0; % diagnostic plots.
onlyPlotMean=0; % then does not plot rasters, just means.

lt_neural_motifRastMult(SongDat, NeurDat, Params, MotifList_regexp, ...
    predur, postdur, LinearWarp, suppressplots, onlyPlotMean)

% === TO DO, LIMIT TO A SINGLE CHANNEL


%% +++++++++++++++++++++++++++++++++++++++++++++++++
%% [debug] extract data for single channel of single file

filename='bk7_160809_110736.rhd';
chan=2;

[amplifier_data, ~, frequency_parameters, ...
    board_adc_data, board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(filename);

data=amplifier_data(2,:);

save('data.mat', 'data');

% ===== after using wave_clus
times_data=load('times_data.mat');
numclust=size(times_data.cluster_class,2);

% - plot this file
    ChansToPlot.DigChans_zero=[]; % make string "all" to plot all that exist. empty array to ignore
    ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
    ChansToPlot.AmpChans_zero=amplifier_channels(chan).chip_channel;
    % ChansToPlot.AmpChans_zero=[10 14 18 23];

    neuralFiltLow=500;

    PlotWhat.raw=0;
    PlotWhat.filt=1;
    PlotWhat.rect_sm=1;
    PlotWhat.raster=1;

    Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
    Raster.ThrXNoise=times_data.par.stdmin; % threshold for units, used for all channels, uses absolute data for peak detection
    Raster.PosOrNeg=-1; % -1 or +1, for crossings.
lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster)

% - overlay spikes detected from wave_clus
plotcols={'g', 'b', 'm', 'c'};
for i=1:numclust
   inds=find([times_data.cluster_class(:,1)]==i);
   
   for ii=inds
      
       spktime=times_data.cluster_class(ii, 2); % in ms
       
       line([spktime spktime], [-20 -40], 'Color', plotcols{i}, 'LineWidth',2);
   end
end




%% output
% for each song file, have a data structure containing various things







