function lt_neural_v2_LEARNING_MetaSummary

%% PARAMs


%% load summary struct

BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 1;
BatchesDesired = {};
ChannelsDesired = [];
[~, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);


%% 
% ========= EXTRACT LEARNING EXPERIMENT INFORMATION
LearnMetaDat = lt_neural_v2_LoadLearnMetadat;

NumBirds = length(LearnMetaDat.bird);

for i=1:NumBirds
   
    birdname = LearnMetaDat.bird(i).birdname;
    
     ListOfExpts = LearnMetaDat.bird(i).info(1,:);
     ListOfExpts = unique(ListOfExpts(~isempty(ListOfExpts)));
    
    for ii=1:length(ListOfExpts)
       
        exptname = ListOfExpts{ii};
        inds = strcmp(LearnMetaDat.bird(i).info(1,:), exptname);
        TargSyls = LearnMetaDat.bird(i).info(2,inds);
        Transitions = LearnMetaDat.bird(i).info(3:end, inds);
        
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         title([birdname '-' exptname]);
        
        % =========== FIND DATA FOR THIS EXPERIMENT
        birdind = strcmp({SummaryStruct.birds.birdname}, birdname);
        assert(sum(birdind)==1, 'weird');
        
        exptind = strcmp({SummaryStruct.birds(birdind).neurons.exptID}, exptname);
        
        if ~any(exptind)
           
            % plot blank - no data for some reason
            lt_plot_text(0, 0.5, 'no neurons found', 'r');
             
            continue
        end
        
        summarystruct_tmp = struct;
        summarystruct_tmp.birds(1).birdname = birdname;
        summarystruct_tmp.birds(1).neurons = SummaryStruct.birds(birdind).neurons(exptind);

        collectWNhit = 0;
        onlyCollectTargSyl = 1;
        MOTIFSTATS = lt_neural_v2_ANALY_ExtractMotif(summarystruct_tmp, collectWNhit, onlyCollectTargSyl);
        
        % --- enter targ syls
        tmp = lt_neural_v2_LoadLearnMetadat;
        indbird = strcmp({tmp.bird.birdname}, birdname);
        indexpt = strcmp(tmp.bird(indbird).info(1,:), exptname);
        TargSyls = tmp.bird(indbird).info(2,indexpt);
        MOTIFSTATS.params.TargSyls = TargSyls;
        
        lt_neural_v2_ANALY_LearningPlot1(summarystruct_tmp, MOTIFSTATS);
        title([birdname '-' exptname '-' TargSyls{:}])
    end
    
end

autoArrangeFigures
