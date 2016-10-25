%% lt 8/19/16 - % EXTRACT SONG, LABEL, ONSETS, SPIKE DATA
% Run from day folder. run after performing wave_clus and verifying
% results.

function [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board)

% channel_board=14; neural
% batchf='BatchTest'; the file used when first concat neural data for
% the above must be the names used for original concatenation - i.e.
% wave_clus folder name

% -- extract labels
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)

%% EXTRACT AND CONCAT ALL SONG FILES and neural files
load('MetaDat.mat'); % contains file names
% --- load concat neural and spikes
% neural_cat=load('data.mat');
spikes_cat=load('times_data.mat');
% fs=metaDat(1).fs; % assumes all have same fs, as they should.

% --- load audio and old neural (from actual files)
cd ..

AllSongs=[];
% AllNeural_old=[];
AllLabels=[];
AllOnsets=[];
AllOffsets=[];
cumulative_filedur=0; % keeps track as concatenating

for i=1:length(metaDat)
    
    if isfield(metaDat, 'songDat')
        % then don't need to reload raw
        AllSongs=[AllSongs metaDat(i).songDat];
    else
        % relaod raw
        
        % -- load original sound and neural
        [amplifier_data,~,frequency_parameters, board_adc_data, ...
            board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(i).filename);
%         ind=find([amplifier_channels.chip_channel]==channel_board);
        
        AllSongs=[AllSongs board_adc_data(1,:)];
        
        %    AllNeural_old=[AllNeural_old amplifier_data(ind, :)];
    end
    
    % -- load labels, onsets, offsets
    tmp=load([metaDat(i).filename '.not.mat']);
    AllLabels=[AllLabels tmp.labels];
    
    % convert onsets to onset relative to start of entire concatenated file
    onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
    offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
    
    AllOnsets=[AllOnsets onsets_cum'];
    AllOffsets=[AllOffsets offsets_cum'];
    
    
    % duration of this song file (sec)
    filedur=metaDat(i).numSamps/metaDat(i).fs;
    cumulative_filedur=cumulative_filedur + filedur;
end



SongDat.AllSongs=AllSongs;
SongDat.AllLabels=AllLabels;
SongDat.AllOnsets=AllOnsets;
SongDat.AllOffsets=AllOffsets;
NeurDat.spikes_cat=spikes_cat;
NeurDat.metaDat=metaDat;
Params.batchf=batchf;
Params.channel_board=channel_board;


% ---- sanity check, plot onsets, offsets, cum times
% figcount=1;
% subplotrows=4;
% subplotcols=1;
% fignums_alreadyused=[];
% hfigs=[];
% tt=[1:length(AllSongs_old)]/frequency_parameters.amplifier_sample_rate;
% hsplots=[];
% 
% % - a. song raw
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% plot(tt, AllSongs_old);
% hsplots=[hsplots hsplot];
% 
% % - b. song spec
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% lt_plot_spectrogram(AllSongs_old, metaDat(1).fs, '', 0)
% hsplots=[hsplots hsplot];
% 
% % onsets and labels
% for i=1:length(AllOnsets)
%     line([AllOnsets(i) AllOffsets(i)], [0 0], 'LineWidth', 2);
%     lt_plot_text(AllOnsets(i), 100, AllLabels(i), 'r')
% end

disp('DONE ! ...');

