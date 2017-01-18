function [ANOVA_STRUCT] = lt_neural_MultNeur_AnovaTcourse_sub1(GroupLabels, ...
    GroupSampleSizes, predur, postdur, SongDat, NeurDat, NeuronDatabase, nn, ...
    min_trials_for_FR_and_ANOVA)



%% 
for ii=1:length(GroupLabels)
    
    regexpr_str = GroupLabels{ii};
%     predur = 0.15;
%     postdur = 0.02;
    alignByOnset = 1;
    
    % - Extract data
    [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, '', ...
        regexpr_str, predur, postdur, alignByOnset, '');
    
    % ---- TAKE RANDOM SUBSET IF NEEDED
    if ~isnan(GroupSampleSizes(ii))
        % then take random subset equal to entry
        % IN PROGRESS
        
        inds = randperm(length(SegmentsExtract));
        inds = sort(inds(1:GroupSampleSizes(ii)));
        SegmentsExtract = SegmentsExtract(inds);
        
    end
    
    
    % ---- throw out if any transitions too long
    % BELOW COPIED, CHANGE TO FIT THIS CODE
%             tokeninds=[SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis];
%             
%             if strcmp(branchtype, 'conv');
%                 tmp = SongDat.AllOnsets(tokeninds) - SongDat.AllOffsets(tokeninds-1);
%                 inds_to_remove = tmp > maxTransTime;
%             elseif strcmp(branchtype, 'div')
%                 tmp = SongDat.AllOnsets(tokeninds+1) - SongDat.AllOffsets(tokeninds);
%                 inds_to_remove = tmp > maxTransTime;
%             end
%             if sum(inds_to_remove) > 1
%                 disp(['Removed ' num2str(sum(inds_to_remove)) ' transtiions (too long)']);
%             end
%             SegmentsExtract(inds_to_remove) = [];
             
    % ------ GET TIMING OF PRE AND POST SYL - COPIED CODE, CHANGE TO FIT
    % THIS CODE
%                 % ---- GET TIMING OF PRECEDING SYL
%             globalind_presyl = [SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]-1;
%             onsets_presyl = SongDat.AllOnsets(globalind_presyl);
%             offsets_presyl = SongDat.AllOffsets(globalind_presyl);
%             
%             onsets_presyl = median(onsets_presyl - ...
%                 [SegmentsExtract.global_ontime_motifInclFlank]); % get relative to onset of extracted data
%             offsets_presyl = median(offsets_presyl - ...
%                 [SegmentsExtract.global_ontime_motifInclFlank]);
%             
%             
%             % ---- GET TIMING OF FOLLOWING SYL
%             globalind_postsyl = [SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]+1;
%             globalind_postsyl(globalind_postsyl>length(SongDat.AllOnsets))=[]; % remove things that don't have postsyl
%             onsets_postsyl = SongDat.AllOnsets(globalind_postsyl);
%             offsets_postsyl = SongDat.AllOffsets(globalind_postsyl);
%             
%             if isempty(onsets_postsyl)
%                 onsets_postsyl = nan;
%                 offsets_postsyl = nan;
%             else
%             onsets_postsyl = median(onsets_postsyl - ...
%                 [SegmentsExtract(1:length(onsets_postsyl)).global_ontime_motifInclFlank]);
%             offsets_postsyl = median(offsets_postsyl - ...
%                 [SegmentsExtract(1:length(offsets_postsyl)).global_ontime_motifInclFlank]);
%             end
% 
    
    % ---- keep only spikes of a certain cluster number
    for iii = 1:length(SegmentsExtract)
        inds = SegmentsExtract(iii).spk_Clust == NeuronDatabase.neurons(nn).clustnum;
        SegmentsExtract(iii).spk_Clust = SegmentsExtract(iii).spk_Clust(inds);
        SegmentsExtract(iii).spk_Times = SegmentsExtract(iii).spk_Times(inds);
    end
    
    
    % -------------- KEEP ONLY IF HAVE ENOUGH DATA
    if length(SegmentsExtract)<min_trials_for_FR_and_ANOVA
        continue
    end
    
    
    % ----- COLLECT ALL (for this master) INTO ONE STRUCTURE
    % important information for performing ANOVA - i) alignment
    % time, ii) category (for each trial)
    tmp_struct = SegmentsExtract; % only keep the essentials here
    tmp_struct = rmfield(tmp_struct, {'song_filename', 'song_datenum', ...
        'song_ind_in_batch', 'global_ontime_motifInclFlank', 'global_offtime_motifInclFlank', ...
        'fs', 'global_startind_motifonly', 'global_endind_motifonly', ...
        'global_tokenind_DatAlignedToOnsetOfThis', 'spk_Clust'});
    
    if ~exist('ANOVA_STRUCT', 'var')
        ANOVA_STRUCT = tmp_struct;
    else
        ANOVA_STRUCT = [ANOVA_STRUCT tmp_struct];
    end
    
end
