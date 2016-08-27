%% lt 8/18/16 - input regexpt, outputs segments,each of which is the regexp +/- duration.

function [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, regexpr_str, predur, postdur)


% regexpr_str='ghh';
% predur=4; % sec
% postdur=4; % sec
% SongDat, NeurDat, Params, shoudl correspond to same dataset, extracted using lt_neural_ExtractDat

% note: aligns all by the first match ind (byt he onset)



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

%% ------ find motifs
% MODIFIED TO ALIGN TO ONSET OF TOKEN
[startinds, endinds, matchlabs, tokenExtents]=regexp(AllLabels, regexpr_str, 'start', 'end', ...
    'match', 'tokenExtents');

SegmentsExtract=struct;

% -- for each match ind, extract audio + spikes
for i=1:length(tokenExtents)
    
    % on time
    ind=tokenExtents{i}(1);    
    ontime=AllOnsets(ind); % sec
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
    SegmentsExtract(i).global_tokenind_DatAlignedToOnsetOfThis=tokenExtents{i};
end

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
