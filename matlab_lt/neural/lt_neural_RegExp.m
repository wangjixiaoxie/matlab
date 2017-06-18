%% lt 8/18/16 - input regexpt, outputs segments,each of which is the regexp +/- duration.

function [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
    regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
    keepRawSongDat, suppressout, collectWNhit, collectWholeBoutPosition)
%% note on FF stuff
% WNhit collection is only performed if FF collection is performed - can
% modify easily to make them independent. Needs raw audio to extract WNhit
% data...

%% note 3/21/17
% Might miss some whole bouts. reason is potential motifs are first detected by
% the onset_pre_threshold criterion (which is about 1s). then I ask how many of
% those are preceded and followed by gap larger than WHOLEBOUTS_edgedur, which
% is larger than onset_pre_threshold. So some will fail. If want all that pass
% onset_pre_threshold to be included, then make  onset_pre_threshold = WHOLEBOUTS_edgedur;
% (i.e. make WHOLEBOUTS_edgedur shorter). problem with that is that then might
% get some bouts with some calls ppreceding and following motif.

%% keepRawSongDat

% will keep raw dat if ask for FF.
% will discard if don't ask, (unless set keepRawSongDat = 1 [default = 0]);

if ~exist('keepRawSongDat', 'var')
    keepRawSongDat = 0;
end

if ~exist('suppressout', 'var')
    suppressout = 1;
end

if ~exist('collectWNhit', 'var')
    collectWNhit=1;
end

if ~exist('collectWholeBoutPosition', 'var')
    collectWholeBoutPosition=0; % to get position of a given datapoint within its bout
end

%% --
% FFparams.collectFF=1;
% FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
%             'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
% FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
%             'v', [2450 4300]};
% FFparams.FF_PosRelToken=1; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
%     % +1 is 1 after token
% FFparams.FF_sylName='b'; % Optional: what syl do you expect this to be? if incompatible will raise error
%     % not required (can set as []);



%%
% regexpr_str='ghh'; % if 'WHOLEBOUTS' then will automatically extract
% bouts by onset and offsets.
% predur=4; % sec, time before token onset (or offset, depending on
% alignByOnset)
% postdur=4; % sec
% SongDat, NeurDat, Params, shoudl correspond to same dataset, extracted using lt_neural_ExtractDat
% alignByOnset=1, then aligns all by the onset of token ind. if 0, then by offset of token ind.

%     WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
%     % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS


if ~exist('FFparams', 'var')
    FFparams.collectFF=0;
end

if isempty(FFparams)
    FFparams.collectFF=0;
end


if ~exist('WHOLEBOUTS_edgedur', 'var');
    WHOLEBOUTS_edgedur=[];
end

if ~exist('alignByOnset', 'var');
    alignByOnset=1;
end

UseLastSylAsToken=0;

%% === extract stuff from inputs


AllLabels=SongDat.AllLabels;
AllOnsets=SongDat.AllOnsets;
AllOffsets=SongDat.AllOffsets;
% spikes_cat=NeurDat.spikes_cat;
% metaDat=NeurDat.metaDat;

fs=NeurDat.metaDat(1).fs;

Params.REGEXP.predur=predur;
Params.REGEXP.postdur=postdur;

% PARAMS ONLY USED FOR WHOLEBOUTS EXTRACTION
%      onset_pre_threshold=2; % OLD (changed on 3/21/17, based on wh6pk36)
onset_pre_threshold=1; % in seconds, minimum time preceding and following bout [DO NOT CHANGE! THIS IS FOR "DEFINING" SONG BOUT]
% If you only want to keep those with a long enough quiet time, then do
% following

%     offset_post_threshold=3; % in seconds
min_syls_in_bout=10;


TotalSamps = sum([NeurDat.metaDat.numSamps]);
TotalDurSec = TotalSamps/fs;



%% ------ find motifs

SegmentsExtract=struct;

if strcmp(regexpr_str, 'WHOLEBOUTS')
    gapdurs=[AllOnsets(1) AllOnsets(2:end)-AllOffsets(1:end-1) TotalDurSec-AllOffsets(end)]; % first gap is just dur from start of file to first onset
    
    % - a motif is
    potentialOnsets=find(gapdurs>onset_pre_threshold); % ind of onset syl
    
    inds_tmp=find(diff(potentialOnsets)>=min_syls_in_bout); % indexes thos bouts that pass min syls criterion
    
    bout_firstsyls=potentialOnsets(inds_tmp);
    bout_lastsyls=potentialOnsets(inds_tmp+1)-1; % minus 1 because want to end of the current bout, not the start of the next bout.
    
    bout_firstsyls_beforeflankfilter = bout_firstsyls;
    bout_lastsyls_beforeflankfilter = bout_lastsyls;
    
    
    % ------------
    % === ONLY KEEP THOSE WHICH HAVE FLANKING GAP DURS LONGER THAN
    % CRITERION
    if ~isempty(WHOLEBOUTS_edgedur)
        gapdurs_subset=gapdurs(gapdurs>onset_pre_threshold);
        inds_tmp2=gapdurs_subset(inds_tmp)>WHOLEBOUTS_edgedur ...
            & gapdurs_subset(inds_tmp+1)>WHOLEBOUTS_edgedur;
        
        bout_firstsyls=bout_firstsyls(inds_tmp2);
        bout_lastsyls=bout_lastsyls(inds_tmp2);
    else
        disp('PROBLEM - WHOLEBOUTS_edgedur is empty!!!');
    end
    
    % ============================== TROUBLESHOOTING - PLOT ALL ONSETS IN A
    % LINE PLOT
    if (0)
        % 1) song
        hsplots = [];
        lt_figure; hold on;
        
        hsplot= lt_subplot(2,1,1); hold on;
        plot([1:length(SongDat.AllSongs)]./NeurDat.metaDat(1).fs, ...
            SongDat.AllSongs, 'k');
        title('song');
        hsplots = [hsplots hsplot];
        
        % 2) syls
        hsplot = lt_subplot(2,1,2); hold on;
        for i=1:length(AllOnsets)
            line([AllOnsets(i) AllOffsets(i)], [0 0], 'LineWidth', 2);
        end
        
        for i=1:length(bout_firstsyls)
            line([AllOnsets(bout_firstsyls(i)) AllOnsets(bout_firstsyls(i))]-WHOLEBOUTS_edgedur,...
                ylim, 'Color', 'g');
            line([AllOffsets(bout_lastsyls(i)) AllOffsets(bout_lastsyls(i))]+WHOLEBOUTS_edgedur, ...
                ylim , 'Color', 'r');
        end
        
        for i=1:length(bout_firstsyls_beforeflankfilter)
            line([AllOnsets(bout_firstsyls_beforeflankfilter(i)) AllOnsets(bout_firstsyls_beforeflankfilter(i))],...
                ylim, 'Color', 'b');
            line([AllOffsets(bout_lastsyls_beforeflankfilter(i)) AllOffsets(bout_lastsyls_beforeflankfilter(i))], ...
                ylim , 'Color', 'b');
        end
        hsplots = [hsplots hsplot];
        
        title('blue: potential motifs...redgreen: those that pass edge_dur test');
        xlabel('time of syl (s)');
        
        linkaxes(hsplots, 'x');
    end
    
    % ----------------------
    
    
    % ---- convert to variables used below
    tokenExtents=bout_firstsyls;
    startinds=bout_firstsyls;
    endinds=bout_lastsyls;
    matchlabs={};
    for j=1:length(bout_firstsyls)
        matchlabs{j}=AllLabels(bout_firstsyls(j):bout_lastsyls(j));
    end
    
    % --- if want to use last syl as token.
    if UseLastSylAsToken==1;
        tokenExtents=bout_lastsyls;
        startinds=bout_lastsyls;
        endinds=bout_lastsyls;
        matchlabs={};
        for j=1:length(bout_lastsyls)
            matchlabs{j}=AllLabels(bout_lastsyls(j):bout_lastsyls(j));
        end
    end
else
    % then is a regexp string
    % MODIFIED TO ALIGN TO ONSET OF TOKEN
    [startinds, endinds, matchlabs, tokenExtents]=regexp(AllLabels, regexpr_str, 'start', 'end', ...
        'match', 'tokenExtents');
    
    functmp = @(X) X(1);
    tokenExtents=cellfun(functmp, tokenExtents); % convert from cell to vector.
end


%% =================== IF WANT TO GET WHOLEBOUTS INFO TO GET DATAPOINT
% POSITION IN BOUT (BUT NOT ACTUALLY EXTRACT WHOLEBOUTS AS MAIN DATAPOINT)
if collectWholeBoutPosition==1 & strcmp(regexpr_str, 'WHOLEBOUTS') ==0
    gapdurs=[AllOnsets(1) AllOnsets(2:end)-AllOffsets(1:end-1) TotalDurSec-AllOffsets(end)]; % first gap is just dur from start of file to first onset; % last gap is just offset to end of file
    
    % - a motif is
    potentialOnsets=find(gapdurs>onset_pre_threshold); % ind of onset syl
    
    inds_tmp=find(diff(potentialOnsets)>=min_syls_in_bout); % indexes thos bouts that pass min syls criterion
    
    wholebout_firstsyls=potentialOnsets(inds_tmp);
    wholebout_lastsyls=potentialOnsets(inds_tmp+1)-1; % minus 1 because want to end of the current bout, not the start of the next bout.
    
    % ---- get match labels
    wholebout_matchlabs={};
    for j=1:length(wholebout_firstsyls)
        wholebout_matchlabs{j}=AllLabels(wholebout_firstsyls(j):wholebout_lastsyls(j));
    end
    
    % TROUBLESHOOT
    if (0)
        % 1) song
        hsplots = [];
        lt_figure; hold on;
        
        %         hsplot= lt_subplot(2,1,1); hold on;
        %         plot([1:length(SongDat.AllSongs)]./NeurDat.metaDat(1).fs, ...
        %             SongDat.AllSongs, 'k');
        %         title('song');
        %         hsplots = [hsplots hsplot];
        
        % 2) syls
        hsplot = lt_subplot(2,1,2); hold on;
        for i=1:length(AllOnsets)
            line([AllOnsets(i) AllOffsets(i)], [0 0], 'LineWidth', 2);
            %             lt_plot_text(AllOnsets(i), 0.01, AllLabels(i));
        end
        
        %         for i=tokenExtents
        %             lt_plot_text(AllOnsets(i), 0.01, AllLabels(i));
        %         end
        for i=1:length(wholebout_firstsyls)
            line([AllOnsets(wholebout_firstsyls(i)) AllOnsets(wholebout_firstsyls(i))],...
                ylim, 'Color', 'g');
            line([AllOffsets(wholebout_lastsyls(i)) AllOffsets(wholebout_lastsyls(i))], ...
                ylim , 'Color', 'r');
        end
        
        xlabel('time of syl (s)');
        
        plot(AllOnsets(638), 0, 'ok')
    end
    
end

%%
HitSyls_TEMP={};
MisSyls_TEMP={};

BoutNumsAll = [];
PosInBoutAll = [];
RendInBoutAll = [];

% -- for each match ind, extract audio + spikes
for i=1:length(tokenExtents)
    
    % on time
    ind=tokenExtents(i);
    if alignByOnset==1
        ontime=AllOnsets(ind); % sec
    else
        % align by offset of token syl
        ontime=AllOffsets(ind); % sec
    end
    
    
    ontime=ontime-predur; % - adjust on time
    onsamp=round(ontime*fs);
    
    % off time
    ind=endinds(i);
    offtime=AllOffsets(ind);
    offtime=offtime+postdur;
    offsamp=round(offtime*fs);
    
    
    spkinds=(NeurDat.spikes_cat.cluster_class(:,2) > ontime*1000) & ...
        (NeurDat.spikes_cat.cluster_class(:,2) < offtime*1000);
    spk_ClustTimes = NeurDat.spikes_cat.cluster_class(spkinds, :); % in sec, relative to onset of the segment
    
    
    if keepRawSongDat ==1
        assert(isfield(SongDat, 'AllSongs'), 'PROBLEM - need to extract songdat before running this');
        
        % this effectively does nothing if also collecting FF.
%         AllSongs=SongDat.AllSongs;
        songseg=SongDat.AllSongs(onsamp:offsamp);
        SegmentsExtract(i).songdat=songseg;
    end
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Extract FF of a specific syllable [OPTIONAL]
    if FFparams.collectFF==1
        
        PC=nan;
        FF=nan;
        T=nan;
        
        if isfield(SongDat, 'FFvals');
            %         if (1)
            % -- take old vals
            FF = SongDat.FFvals(tokenExtents(i));
            
        else
            %% === OLD METHOD (CALCULATE FROM RAW AUDIO HERE)
            disp('HAVE TO DO FFVALS MANUALLY :( ');
            % - collect
%             AllSongs=SongDat.AllSongs;
            songseg=SongDat.AllSongs(onsamp:offsamp);
            
            FF_PosRelToken=FFparams.FF_PosRelToken;
            FF_sylName=FFparams.FF_sylName;
            cell_of_freqwinds=FFparams.cell_of_freqwinds;
            cell_of_FFtimebins=FFparams.cell_of_FFtimebins;
            
            prepad_forFF=0.015; postpad_forFF=0.015; % don't change, as this affects temporal window for FF
            
            indForFF=tokenExtents(i)+FF_PosRelToken;
            sylForFF=AllLabels(indForFF);
            
            ontime_forFF=AllOnsets(indForFF);
            ontime_forFF=ontime_forFF-prepad_forFF;
            onsamp_forFF=round(ontime_forFF*fs);
            
            offtime_forFF=AllOffsets(indForFF);
            offtime_forFF=offtime_forFF+postpad_forFF;
            offsamp_forFF=round(offtime_forFF*fs);
            
            songseg_forFF=SongDat.AllSongs(onsamp_forFF:offsamp_forFF);
            
            % -- debug
            if (0)
                lt_figure; hold on;
                lt_plot_spectrogram(songseg_forFF, fs, 1,0);
                pause
                close all;
            end
            
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXTRACT FF
            % -- defaults, sets output (e.g., ff) to this if don't collect
            % actual FF
            collectFF=1;
            
            % ---------------- GO THROUGH POTENTIAL REASONS TO NOT COLLECT FF
            if ~isempty(FF_sylName)
                %             assert(strcmp(FF_sylName, sylForFF)==1, 'Syl is not the desired syl!! - stopping');
                if strcmp(FF_sylName, sylForFF)==0
                    collectFF=0;
                end
            end
            
            indtmp=find(strcmp(sylForFF, cell_of_freqwinds));
            if isempty(indtmp)
                % tel;l use
                disp('EMPTY FREQ WINDOW!! - giving nan for FF');
                collectFF=0;
            elseif isempty(cell_of_freqwinds{indtmp+1})
                disp('EMPTY TIME WINDOW!! - giving nan for FF');
                collectFF=0;
            end
            
            % ------ COLLECT
            if collectFF==1
                indtmp=find(strcmp(sylForFF, cell_of_freqwinds));
                F_high=cell_of_freqwinds{indtmp+1}(2);
                F_low=cell_of_freqwinds{indtmp+1}(1);
                
                indtmp=find(strcmp(cell_of_FFtimebins, sylForFF));
                mintime=cell_of_FFtimebins{indtmp+1}(1); % sec
                maxtime=cell_of_FFtimebins{indtmp+1}(2);
                
                [FF, PC, T]= lt_calc_FF(songseg_forFF, fs, [F_low F_high], [mintime maxtime]);
            end
        end
        
        SegmentsExtract(i).FF_val=FF;
        %% ===================================================
        % FIGURE OUT IF WN HIT ON THIS TRIAL (based on clipping of sound)
        % - collect
        if collectWNhit==1
%             AllSongs=SongDat.AllSongs;
            songseg=SongDat.AllSongs(onsamp:offsamp);
            
            FF_PosRelToken=FFparams.FF_PosRelToken;
            prepad_forFF=0.015; postpad_forFF=0.015; % don't change, as this affects temporal window for FF
            
            indForFF=tokenExtents(i)+FF_PosRelToken; % can use this to look for WN in flanking syls.
            
            ontime_forFF=AllOnsets(indForFF);
            ontime_forFF=ontime_forFF-prepad_forFF;
            onsamp_forFF=round(ontime_forFF*fs);
            
            offtime_forFF=AllOffsets(indForFF);
            offtime_forFF=offtime_forFF+postpad_forFF;
            offsamp_forFF=round(offtime_forFF*fs);
            
            songseg_forFF=SongDat.AllSongs(onsamp_forFF:offsamp_forFF);
            
            
            wasTrialHit=[];
            WNonset=[];
            WNoffset=[];
            if max(songseg_forFF)>3.2
                % --- debug, plot spectrogram and sound file for all hits
                %             figure; hold on;
                %             lt_subplot(1,2,1); hold on; plot(songseg_forFF);
                %             lt_subplot(1,2,2); hold on; lt_plot_spectrogram(songseg_forFF, fs, 1, 0);
                %             pause
                %             close all;
                % ----
                wasTrialHit=1;
                WNonset=find(songseg>3.2, 1, 'first'); % timepoint of hit (for entire segment)
                WNoffset=find(songseg>3.2, 1, 'last');
                
                HitSyls_TEMP=[HitSyls_TEMP songseg_forFF];
            else
                wasTrialHit=0;
                MisSyls_TEMP=[MisSyls_TEMP songseg_forFF];
            end
            
            SegmentsExtract(i).hit_WN=wasTrialHit;
            SegmentsExtract(i).WNonset_sec=WNonset/fs;
            SegmentsExtract(i).WNoffset_sec=WNoffset/fs;
            
            
            %         SegmentsExtract(i).FF_pitchcontour=PC;
            %         SegmentsExtract(i).FF_timebase=T;
        end
    end
    
    % =============== FIGURE OUT POSITION OF MOTIF WITHIN ITS BOUT
    if collectWholeBoutPosition==1
        ind; % current syl posotion
        
        boutnum = find(wholebout_firstsyls<=ind & wholebout_lastsyls>=ind);
        if length(boutnum)==0
            SegmentsExtract(i).BOUT_boutnum = nan;
            SegmentsExtract(i).BOUT_posOfTokenInBout = nan;
            SegmentsExtract(i).BOUT_RendInBout = nan;
            
        else
            % assert(length(boutnum)==1, 'sdafasd');
            positionOfTokenInBout = ind - wholebout_firstsyls(boutnum) + 1;
            RendInBout = sum(BoutNumsAll == boutnum)+1;
            
            BoutNumsAll = [BoutNumsAll boutnum];
            PosInBoutAll = [PosInBoutAll positionOfTokenInBout];
            RendInBoutAll = [RendInBoutAll RendInBout];
            
            SegmentsExtract(i).BOUT_boutnum = boutnum;
            SegmentsExtract(i).BOUT_posOfTokenInBout = positionOfTokenInBout;
            SegmentsExtract(i).BOUT_RendInBout = RendInBout;
        end
    end
    
    
    % =================== GET ONSETS/OFFSETS OF ALL SYLS IN MOTIF, RELATIVE
    % TO TIMING OF EACH MOTIF.
    if alignByOnset==1
        %     OnsetsRelTokenOnset = AllOnsets - AllOnsets(tokenExtents(i))
        
        % get only those that are within bounds of this segment's data
        inds = AllOnsets>ontime & AllOnsets<offtime;
        ontimesWithinData = AllOnsets(inds);
        ontimesRelStartDat = ontimesWithinData - ontime; % relative to start of segment
        
        inds = AllOffsets>ontime & AllOffsets<offtime;
        offtimesWithinData = AllOffsets(inds);
        offtimesWithinData = offtimesWithinData - ontime; % relative to start of segment
        
        % --- compare onsets to acoustic data
        if (0)
            lt_figure; hold on;
            lt_plot(ontimesRelStartDat, 0, {'Color', 'g'});
            lt_plot(offtimesWithinData, 0, {'Color', 'r'});
            plot([1:length(songseg)]./fs, (songseg-mean(songseg)).^2, 'k');
            line([predur predur], ylim, 'Color', 'b');
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ========= TIME OF SONG FILE FOR THIS SEGMENT
    globalOnsetTime=AllOnsets(tokenExtents(i)); % sec
    globalOnsetSamp=globalOnsetTime*fs;
    cumOnsSamps=cumsum([0 NeurDat.metaDat.numSamps]);
    songind=find((globalOnsetSamp-cumOnsSamps)>0, 1, 'last');
    %     disp(['songind = ' num2str(songind)]);
    songfname=NeurDat.metaDat(songind).filename;
    
    SegmentsExtract(i).song_filename=songfname;
    SegmentsExtract(i).song_datenum=lt_neural_fn2datenum(songfname);
    SegmentsExtract(i).song_ind_in_batch=songind;
    
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    SegmentsExtract(i).spk_Clust=spk_ClustTimes(:,1)';
    SegmentsExtract(i).spk_Times=(spk_ClustTimes(:,2)/1000)'-ontime;
    SegmentsExtract(i).global_ontime_motifInclFlank=ontime;
    SegmentsExtract(i).global_offtime_motifInclFlank=offtime;
    SegmentsExtract(i).matchlabel=matchlabs{i};
    SegmentsExtract(i).fs=fs;
    SegmentsExtract(i).global_startind_motifonly=startinds(i);
    SegmentsExtract(i).global_endind_motifonly=endinds(i);
    SegmentsExtract(i).global_tokenind_DatAlignedToOnsetOfThis=tokenExtents(i);
    
    actualmotifdur=(offtime - postdur) - (ontime + predur);
    SegmentsExtract(i).actualmotifdur=actualmotifdur;
    
    SegmentsExtract(i).sylOnTimes_RelDataOnset=ontimesRelStartDat;
    SegmentsExtract(i).sylOffTimes_RelDataOnset=offtimesWithinData;
    
end

if suppressout==0
    disp(['DONE! EXTRACTED ' num2str(length(SegmentsExtract)) ' segments']);
end



%% ==== DEBUG - CHECK HIT DETECTION
if suppressout==0
    if ~isempty(HitSyls_TEMP) | ~isempty(MisSyls_TEMP)
        lt_figure; hold on;
        lt_subplot(1,2,1); hold on; title('hits');
        for j=1:length(HitSyls_TEMP)
            plot(HitSyls_TEMP{j}, ':', 'Color', [rand rand rand]);
        end
        lt_subplot(1,2,2); hold on; title('misses');
        for j=1:length(MisSyls_TEMP)
            plot(MisSyls_TEMP{j}, ':', 'Color', [rand rand rand]);
        end
        
        lt_subtitle(regexpr_str);
    end
end

