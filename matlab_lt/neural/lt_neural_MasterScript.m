
%% ===== plot song

songname='bk7_160809_110736.rhd';
songnotmat='bk7_160809_110736.rhd.not.mat';

lt_neural_extractaudio(songname);

% ---- overlay onsets and offsets [first label with evsonganaly]
notmatstruct=load(songnotmat);

for i=1:length(notmatstruct.onsets)
    
    line([notmatstruct.onsets(i) notmatstruct.offsets(i)]./1000, [1.65 1.65], 'Color','r','LineWidth', 7);
end


%% get a list of sampling rates of all files


%% ====== plot single file dat [align neural and song]
close all;
filename='bk7_160809_112713.rhd';
ChansToPlot.DigChans_zero=[]; % make string "all" to plot all that exist. empty array to ignore
ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
ChansToPlot.AmpChans_zero=[9 14 19];
% ChansToPlot.AmpChans_zero=[10 14 18 23];

neuralFiltLow=500;

PlotWhat.raw=0;
PlotWhat.filt=1;
PlotWhat.rect_sm=1;
PlotWhat.raster=1;

Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
Raster.ThrXNoise=6; % threshold for units, used for all channels, uses absolute data for peak detection
Raster.PosOrNeg=-1; % -1 or +1, for crossings.
lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster)



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ANALYSIS PIPELINE
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% label songs using modified evsonganaly
clear all; close all;
% batchf='batch_test';
batchf='BatchTest';

%% ==== exploratory - concat all audio and neural and plot for each neural channel
ChansToPlot=[9 14 19];
lt_neural_concatExplore(batchf, ChansToPlot)


%% ==== concatenate multiple files for given channel, [and saves]
% based on expectation of duration for single unit.
% -- saves into downstream folder

channel_board=14;
lt_neural_concatOneChan(batchf, channel_board)

%% ==== run wave_clus on this concatted data


%% ==== [SANITY CHECK] plot song files in entirety, aligned to extracted spikes
% -- makes multiple plots if too much dat.
close all;
lt_neural_AlgnWavclus(batchf, channel_board);

%% ==== EXTRACT SONG, LABEL, ONSETS, SPIKE DATA

[SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board)

%% ===== Align motifs with neural data
regexpr_str='ghh';
predur=3; % sec
postdur=4; % sec
channel_board=14;
batchf='BatchTest';

close all;
lt_neural_motifRaster(batchf, channel_board, regexpr_str, predur, postdur)


%% =========== PLOT onset and offset firing 
% (single unit, multi unit, smoothed mu,...)






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







