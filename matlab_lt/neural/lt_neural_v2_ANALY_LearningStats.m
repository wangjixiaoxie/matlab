%% lt 5/22/17 - across all birds, expts, neurons, plot learning related statistics (see other code
% for visualization of raw data)
function lt_neural_v2_ANALY_LearningStats(MOTIFSTATS_Compiled)

% for smoothing fr
window_neural = 0.0075; % neural;
windshift_neural = 0.002;



%%
NumBirds = length(MOTIFSTATS_Compiled.birds);



%% GET MEAN FR - go thru all experiments

for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    
    for ii=1:numexpts
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        for j=1:nummotifs
            
            % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
            % INTERVALS FOR NEURONS
            
            % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
            numneurons = length(MotifStats.neurons);
            
            for nn=1:numneurons
                
                
                % ==== For each rendition, get smoothed firing rate
                clustnum = MotifStats.neurons(nn).clustnum;
                assert(clustnum == SummaryStruct.birds(1).neurons(nn).clustnum, 'asdfaf');
                
                MotifStats.neurons(nn).motif(j).SegmentsExtract = ...
                    lt_neural_SmoothFR(MotifStats.neurons(nn).motif(j).SegmentsExtract, ...
                    clustnum);
                
                
            end
        end
        
        % === stick back into main structure
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS = MotifStats;
        
    end
end


%% calculate learning metric (e.g. correlation of smoothed FR)

for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    
    for ii=1:numexpts
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        for j=1:nummotifs
            
            lt_figure; hold on;
            
            % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
            % INTERVALS FOR NEURONS
            
            % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
            numneurons = length(MotifStats.neurons);
            
            for nn=1:numneurons
                
                segmentsextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                
                % === 1) GET BASELINE MEAN SMOOTHED RATE
                
                WNonDnum = datenum(SummaryStruct.birds(1).neurons(nn).LEARN_WNonDatestr, ...
                    'ddmmmyyyy-HHMM');
                baseInds = [segmentsextract.song_datenum] < WNonDnum;
                
                
                
                
                
            end
        end
    end
end

