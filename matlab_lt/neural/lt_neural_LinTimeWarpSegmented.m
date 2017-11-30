function segextract = lt_neural_LinTimeWarpSegmented(segextract, regionstowarp)
%% NOTE!!
% this assumes data aligned to syl onset - not sure if will work if aligned
% to offset - best to empirically test.

%%
% regionstowarp = [3 4]; % [1 2 3 ..], which regions to warp in where index refers to, in order, first syl, first gap, 
% second syl second gap , etc.

%%
% segextract = CLASSES.birds(1).neurons(1).branchnum(1).SEGEXTRACT.classnum(1).SegmentsExtract;
% params = CLASSES.birds(1).neurons(1).branchnum(1).SEGEXTRACT.classnum(1).Params;
% SummaryStruct;
%% ======= extract onset and offset of all syllables in motif

if(0)
    lt_figure; hold on;
    for i=1:length(segextract)
        plot(segextract(i).motifsylOnsets, i,'ok');
        plot(segextract(i).motifsylOffsets, i, 'or');
        plot(segextract(i).spk_Times, i, 'xb');
    end
end

%% ======


% ================== get global medians for all REgions (see below for
% definition of regions
RegiondurationAll = nan(length(segextract), 2*length(segextract(1).motifsylOnsets)-1);
SpkRegionsAll = {};
SpkTimeWithinRegAll = {};

for i=1:length(segextract)
    
    % =============== COLLECT Region DURATION (i.e. duration of [syl gap syl
    % gap ...syl], in order
    regionduration = nan(1, 2*length(segextract(i).motifsylOnsets) -1); % syls and gaps
    allsyldurs =  segextract(i).motifsylOffsets-segextract(i).motifsylOnsets;
    allgapdurs = segextract(i).motifsylOnsets(2:end) - segextract(i).motifsylOffsets(1:end-1);
    
    regionduration(1:2:end) = allsyldurs;
    regionduration(2:2:end) = allgapdurs;
    
    RegiondurationAll(i, :) = regionduration;
    
    % ================== COLLECT REGION ONSETS
    regiononsets_relmotif = nan(1, 2*length(segextract(i).motifsylOnsets) -1);
    regiononsets_relmotif(1:2:end) = segextract(i).motifsylOnsets;
    regiononsets_relmotif(2:2:end) = segextract(i).motifsylOffsets(1:end-1);
    
    
    % ================== for each spike, determine the Region it is in, and the
    % position (i.e. time from onset) within the Region
%     assert(length(unique(segextract(i).spk_Clust))==1, 'mult clusters...');
    spktimes = segextract(i).spk_Times;
    spk_regions = [];
    spk_timeWithinReg = [];
    for j=1:length(spktimes)
        thisregion = find(spktimes(j)>regiononsets_relmotif, 1, 'last');
        thisTimeWithinRegion = spktimes(j)-regiononsets_relmotif(thisregion);
        
        assert(thisTimeWithinRegion<regionduration(thisregion), 'problem: if duration greater, then shoud actually be in next region');
        
        spk_regions = [spk_regions thisregion];
        spk_timeWithinReg = [spk_timeWithinReg thisTimeWithinRegion];
    end
    SpkRegionsAll = [SpkRegionsAll spk_regions];
    SpkTimeWithinRegAll = [SpkTimeWithinRegAll spk_timeWithinReg];
    
    % % ------ onset times of all regions
    % RegionOnsets = [0 cumsum(regionduration)];
    % RegionOnsets(end) = [];
    
end

assert(~any(isnan(RegiondurationAll(:))), 'asfasd');

% ==== get global median region durations
RegiondurationMedian = median(RegiondurationAll,1);


% ============================== for each trial, warp to global median
SpkTimeWithinRegAll_scaled = SpkTimeWithinRegAll;
RegiondurationAll_scaled = RegiondurationAll;

for rr = regionstowarp
    for i=1:length(segextract)
    
        scalefactor = RegiondurationMedian(rr)/RegiondurationAll(i,rr); % scale data by this factor
        
        % -------------------- scale region duration (this will trivilaly
        % equal median
        RegiondurationAll_scaled(i,rr) = RegiondurationAll_scaled(i,rr)*scalefactor;
        
        % --------------------- scale all spike times 
        inds = SpkRegionsAll{i}==rr;
        SpkTimeWithinRegAll_scaled{i}(inds) = ...
            SpkTimeWithinRegAll_scaled{i}(inds).*scalefactor;
    end
end


% ============================== reconstruct spike times and put back into
% segextract
for i=1:length(segextract)
   
    regiondurations = RegiondurationAll_scaled(i,:);
    
    regiononsets = [0 cumsum(regiondurations)]; % - convert to onsets
    regiononsets(end) = [];
    
    % ================================================  extract new
    % onsets/offsets
    sylonsets_new = regiononsets(1:2:end);
    syloffsets_new = [regiononsets(2:2:end) sylonsets_new(end)+regiondurations(end)];
    
    
    % ===================================================================
    % --- convert spiketimes from times (relative region onsets) to times
    % (relative onset of entire motif)
    spk_regions = SpkRegionsAll{i};
    spktimes_withinregion = SpkTimeWithinRegAll_scaled{i};
    spktimes_motif = [];
    for j=1:length(spktimes_withinregion)
        thisregion = spk_regions(j);
        thistime = spktimes_withinregion(j);
                
        thistime_global = regiononsets(thisregion)+thistime;
        
        spktimes_motif = [spktimes_motif thistime_global];        
    end
    
    % #################################################################
    % ############ finally, convert from time relative to first syl to time
    % relative to (tokensylonset - motifpredur);
    tokenpos = segextract(i).global_tokenind_DatAlignedToOnsetOfThis - ...
        segextract(i).global_startind_motifonly+1; 
    motifpredur = segextract(i).motifsylOnsets(tokenpos);
    tokenonset = regiononsets(2*tokenpos-1); % picks out onset time of token
    
    % ================ SPKTIMES 
    spktimes_reltokenminuspredur = (spktimes_motif - tokenonset) + motifpredur;
    
    % ================ SYL ONSETS/OFFSETS
    sylonsets_new = (sylonsets_new - tokenonset) + motifpredur;
    syloffsets_new = (syloffsets_new - tokenonset) + motifpredur;
%     assert(all(diff(syloffsets_new)>0), 'asfasd');
%     assert(all(diff(sylonsets_new)>0), 'asfasd');

    % ################################### OUTPUT TO segextract
    segextract(i).LinTWSeg_spkTimes = spktimes_reltokenminuspredur;
    segextract(i).LinTWSeg_sylonsets= sylonsets_new;
    segextract(i).LinTWSeg_syloffsets = syloffsets_new;
            
end

%% ============== plot, compare pre warp to post warp

if(1)
    lt_figure; hold on;
    hsplots =[];
    
    % --- pre warp
    hsplot = lt_subplot(1,2,1); hold on;
    title('prewarp');
    hsplots = [hsplots hsplot];
    for i=1:length(segextract)
        plot(segextract(i).motifsylOnsets, i,'ok');
        plot(segextract(i).motifsylOffsets, i, 'or');
        plot(segextract(i).spk_Times, i, 'xb');
    end
    
    % --- post warp
    hsplot = lt_subplot(1,2,2); hold on;
    title('post warp');
    hsplots = [hsplots hsplot];
    for i=1:length(segextract)
        plot(segextract(i).LinTWSeg_sylonsets, i,'ok');
        plot(segextract(i).LinTWSeg_syloffsets, i, 'or');
        plot(segextract(i).LinTWSeg_spkTimes, i, 'xb');
    end
    
    linkaxes(hsplots, 'xy');
end


