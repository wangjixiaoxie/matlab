%% lt 8/18/16 - input regexpt, outputs segments,each of which is the regexp +/- duration.

function [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur)


% regexpr_str='ghh'; % if 'WHOLEBOUTS' then will automatically extract
% bouts by onset and offsets.
% predur=4; % sec
% postdur=4; % sec
% SongDat, NeurDat, Params, shoudl correspond to same dataset, extracted using lt_neural_ExtractDat
% alignByOnset=1, then aligns all by the onset of token ind. if 0, then by offset of token ind.

%     WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
%     % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS

if ~exist('WHOLEBOUTS_edgedur', 'var');
    WHOLEBOUTS_edgedur=[];
end

if ~exist('alignByOnset', 'var');
    alignByOnset=1;
end

UseLastSylAsToken=0;

%% === extract stuff from inputs


AllSongs=SongDat.AllSongs;
AllLabels=SongDat.AllLabels;
AllOnsets=SongDat.AllOnsets;
AllOffsets=SongDat.AllOffsets;
spikes_cat=NeurDat.spikes_cat;
% metaDat=NeurDat.metaDat;

fs=NeurDat.metaDat(1).fs;

Params.REGEXP.predur=predur;
Params.REGEXP.postdur=postdur;

if strcmp(regexpr_str, 'WHOLEBOUTS')
    
    onset_pre_threshold=2; % in seconds, minimum time preceding and following bout [DO NOT CHANGE! THIS IS FOR "DEFINING" SONG BOUT]
    % If you only want to keep those with a long enough quiet time, then do
    % following
    
    %     offset_post_threshold=3; % in seconds
    min_syls_in_bout=5;
end



%% ------ find motifs

SegmentsExtract=struct;

if strcmp(regexpr_str, 'WHOLEBOUTS')
    gapdurs=[AllOnsets(1) AllOnsets(2:end)-AllOffsets(1:end-1)]; % first gap is just dur from start of file to first onset
    
    % - a motif is
    potentialOnsets=find(gapdurs>onset_pre_threshold); % ind of onset syl
    
    inds_tmp=find(diff(potentialOnsets)>=min_syls_in_bout); % indexes thos bouts that pass min syls criterion
    
    bout_firstsyls=potentialOnsets(inds_tmp);
    bout_lastsyls=potentialOnsets(inds_tmp+1)-1; % minus 1 because want to end of the current bout, not the start of the next bout.
    
    % ------------
    % === ONLY KEEP THOSE WHICH HAVE FLANKING GAP DURS LONGER THAN
    % CRITERION
    if ~isempty(WHOLEBOUTS_edgedur)
        gapdurs_subset=gapdurs(gapdurs>onset_pre_threshold);
        inds_tmp2=gapdurs_subset(inds_tmp)>WHOLEBOUTS_edgedur ...
            & gapdurs_subset(inds_tmp+1)>WHOLEBOUTS_edgedur;
        
        bout_firstsyls=bout_firstsyls(inds_tmp2);
        bout_lastsyls=bout_lastsyls(inds_tmp2);
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
    
    % - collect
    songseg=AllSongs(onsamp:offsamp);
    
    spkinds=(spikes_cat.cluster_class(:,2) > ontime*1000) & ...
        (spikes_cat.cluster_class(:,2) < offtime*1000);
    spk_ClustTimes = spikes_cat.cluster_class(spkinds, :); % in sec, relative to onset of the segment
    
    
    
    SegmentsExtract(i).songdat=songseg;
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
    
    
    
end

disp(['DONE! EXTRACTED ' num2str(length(SegmentsExtract)) ' segments']);


% OLD, using start ind only
% [startinds, endinds, matchlabs, tokenExtents]=regexp(AllLabels, regexpr_str, 'start', 'end', ...
%     'match', 'tokenExtents');
%
% SegmentsExtract=struct;
%
% % -- for each match ind, extract audio + spikes
% for i=1:length(startinds)
%
%     % on time
%     ind=startinds(i);
%     ontime=AllOnsets(ind); % sec
%     ontime=ontime-predur; % - adjust on time
%     onsamp=round(ontime*fs);
%
%     % off time
%     ind=endinds(i);
%     offtime=AllOffsets(ind);
%     offtime=offtime+postdur;
%     offsamp=round(offtime*fs);
%
%     % - collect
%     songseg=AllSongs(onsamp:offsamp);
%
%     spkinds=(spikes_cat.cluster_class(:,2) > ontime*1000) & ...
%         (spikes_cat.cluster_class(:,2) < offtime*1000);
%     spk_ClustTimes = spikes_cat.cluster_class(spkinds, :); % in sec, relative to onset of the segment
%
%     SegmentsExtract(i).songdat=songseg;
%     SegmentsExtract(i).spk_Clust=spk_ClustTimes(:,1)';
%     SegmentsExtract(i).spk_Times=(spk_ClustTimes(:,2)/1000)'-ontime;
%     SegmentsExtract(i).global_ontime_motifInclFlank=ontime;
%     SegmentsExtract(i).global_offtime_motifInclFlank=offtime;
%     SegmentsExtract(i).matchlabel=matchlabs{i};
%     SegmentsExtract(i).fs=fs;
%     SegmentsExtract(i).global_startind_motifonly=startinds(i);
%     SegmentsExtract(i).global_endind_motifonly=endinds(i);
% end



