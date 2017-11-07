%% lt 8/19/16 - % EXTRACT SONG, LABEL, ONSETS, SPIKE DATA
% Run from day folder. run after performing wave_clus and verifying
% results.

function [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, ...
    extractSound, clustToKeep) 
%% lt 11/4/17 - can ask to only keep one clust (default is to take all clusts)
% do not define, or keep empty, to take all clusts.

if ~exist('clustToKeep', 'var')
    clustToKeep = [];
end
    
%% making it optional to extract sound data

if ~exist('extractSound', 'var')
    extractSound = 0;
end


%%
% channel_board=14; neural
% batchf='BatchTest'; the file used when first concat neural data for
% the above must be the names used for original concatenation - i.e.
% wave_clus folder name

% -- extract labels
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)

%% EXTRACT AND CONCAT ALL SONG FILES and neural files
% tic
load('MetaDat.mat'); % contains file names
% toc
% --- load concat neural and spikes
% neural_cat=load('data.mat');
spikes_cat=load('./times_data.mat');
% fs=metaDat(1).fs; % assumes all have same fs, as they should.

% --- load audio and old neural (from actual files)
cd ..

AllSongs=nan(1, sum([metaDat.numSamps]));
% AllNeural_old=[];
AllLabels=[];
AllOnsets=[];
AllOffsets=[];
AllSongNum = []; 
cumulative_filedur=0; % keeps track as concatenating

if extractSound==1
    % look for saved sound file
    soundfilename = ['SongDat_' batchf '.mat'];
    if exist(soundfilename, 'file') ==2
        load(soundfilename); % to SongCellArray
        assert(length(SongCellArray)==length(metaDat), 'problem!');
    else
        disp('NO SAVED SOUND FILE !');
    end
end

% tic
counter = 1;
for i=1:length(metaDat)
    
    if extractSound==1
        
        if exist('SongCellArray', 'var')==1 % ideally use this.
            dattmp = SongCellArray{i};
        else
            if isfield(metaDat, 'songDat') % old version (obsolete after remove songdart from metadat)
                % then don't need to reload raw
                dattmp = metaDat(i).songDat(1,:);
                
            else
                % relaod raw
                
                % -- load original sound and neural
                [~, board_adc_data] = pj_readIntanNoGui_AudioOnly(metaDat(i).filename);
                
                % NEW VERSION
                dattmp = board_adc_data(1,:);
            end
        end
        
        inds = counter:(counter+length(dattmp)-1);
        AllSongs(inds) = dattmp;
        
        counter = counter+length(dattmp);
    end
    
    
    % -- load labels, onsets, offsets
    if exist([metaDat(i).filename '.not.mat'], 'file')==2
        % then file .notmat exists
        tmp=load([metaDat(i).filename '.not.mat']);
        AllLabels=[AllLabels tmp.labels];
        
        % convert onsets to onset relative to start of entire concatenated file
        onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
        offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
        
        AllOnsets=[AllOnsets onsets_cum'];
        AllOffsets=[AllOffsets offsets_cum'];
        AllSongNum = [AllSongNum ones(1, length(tmp.labels))*i];
    else
%         disp(['NOTE: miossing .not.mat for ' metaDat(i).filename ' (SKIPPING)']);
    end
    % duration of this song file (sec)
    filedur=metaDat(i).numSamps/metaDat(i).fs;
    cumulative_filedur=cumulative_filedur + filedur;
    
end
% toc

   
if extractSound==0
    % --- remove sound from metaDat (if it is present)
    if isfield(metaDat, 'songDat');
        metaDat = rmfield(metaDat, 'songDat');
    end
end

% =============== INLCUDE PREVIOUSLY EXTRACTED STATS (E.G. FF)
cd(datdir)
% ------ 1) FF
if exist('extractFF.mat', 'file')
    FFvals = load('extractFF.mat');
    FFparams = load('extractFF_params');
    
    % sanity check - compare params to current params
    if (0) % old version, now just don't extract if fail
        assert(all(FFparams.Params.allLabels == AllLabels), 'PROBLEM');
    end
    %         assert(all((FFparams.Params.AllOnsets - AllOnsets)<0.002), 'PROBLEM');
    tmp =0;
    if length(FFparams.Params.allLabels) == length(AllLabels)
        if all(FFparams.Params.allLabels == AllLabels);
            SongDat.FFvals = FFvals.FFvals;
            tmp = 1;
        end
    end
    if tmp==0;
        disp('PROBLEM, need to re extyract FF');
    end
end

if extractSound==1
SongDat.AllSongs=single(AllSongs);
end


%% ============ keep only a single cluster?

if ~isempty(clustToKeep)
   
    indstokeep = spikes_cat.cluster_class(:,1) == clustToKeep;
    
    assert(any(indstokeep)==1, 'problem - this clust does not exist. ..');
    
    spikes_cat.spikes = spikes_cat.spikes(indstokeep, :);
    spikes_cat.cluster_class = spikes_cat.cluster_class(indstokeep,:);
    spikes_cat.forced = spikes_cat.forced(indstokeep);
    spikes_cat.inspk = spikes_cat.inspk(indstokeep, :);
%     spikes_cat.ipermut = spikes_cat.ipermut(indstokeep);
   
end


%% OUTPUT

SongDat.AllLabels=AllLabels;
SongDat.AllOnsets=AllOnsets;
SongDat.AllOffsets=AllOffsets;
SongDat.AllSongNum=AllSongNum;
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

%% old version saved [right before implemented extractSound as option (i.e. saveds ound outside of metadat)
if(0)
    if extractSound==1
        
        if isfield(metaDat, 'songDat')
            % then don't need to reload raw
            % OLD VERSION
            %         AllSongs=[AllSongs metaDat(i).songDat];
            
            % NEW VERSION
            dattmp = metaDat(i).songDat(1,:);
            
            % inds = counter:(counter+length(dattmp)-1);
            % %         assert(all(isnan(AllSongs(inds))), 'asdfasdf');
            %         AllSongs(inds) = dattmp;
            % %         assert(length(inds) ==length(metaDat(i).songDat(1,:)), 'asdfasdf');
            %
            %         counter = counter+length(dattmp);
            
            %         disp([inds(1) inds(end) length(metaDat(i).songDat(1,:))]);
            
        else
            % relaod raw
            
            % -- load original sound and neural
            [~, board_adc_data] = pj_readIntanNoGui_AudioOnly(metaDat(i).filename);
            %         ind=find([amplifier_channels.chip_channel]==channel_board);
            
            % OLD VERSION
            %         AllSongs=[AllSongs board_adc_data(1,:)];
            
            % NEW VERSION
            dattmp = board_adc_data(1,:);
            
            %         inds = counter:(counter+length(dattmp)-1);
            % %         assert(all(isnan(AllSongs(inds))), 'asdfasdf');
            %         AllSongs(inds) = dattmp;
            % %         assert(length(inds) ==length(board_adc_data(1,:)), 'asdfasdf');
            %
            %         counter = counter+length(dattmp);
            
            %         disp([inds(1) inds(end) length(board_adc_data(1,:))]);
            
            %    AllNeural_old=[AllNeural_old amplifier_data(ind, :)];
        end
        
        
        inds = counter:(counter+length(dattmp)-1);
        %         assert(all(isnan(AllSongs(inds))), 'asdfasdf');
        AllSongs(inds) = dattmp;
        %         assert(length(inds) ==length(board_adc_data(1,:)), 'asdfasdf');
        
        counter = counter+length(dattmp);
    end

end