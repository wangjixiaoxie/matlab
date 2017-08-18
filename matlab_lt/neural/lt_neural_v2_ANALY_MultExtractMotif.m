function MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase)

% LearnKeepOnlyBase for learning, if 1, then keeps only baseline periods (i.e. before WN)
    
%%

if ~exist('collectWNhit', 'var')
   collectWNhit = 1; 
end

if ~exist('LearnKeepOnlyBase', 'var');
    LearnKeepOnlyBase =0;
end

%% lt 6/8/17 - multiple birds, extracts motifs/segments, etc

NumBirds = length(SummaryStruct.birds);
MOTIFSTATS_Compiled = struct; % heirarchy: birds --> expt --> neurons

for i=1:NumBirds

    ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
    birdname = SummaryStruct.birds(i).birdname;
    
    for ll=1:length(ListOfExpts)
        exptname = ListOfExpts{ll};
        
        inds = strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname);
        
        SummaryStruct_tmp = struct;
        SummaryStruct_tmp.birds(1).neurons = SummaryStruct.birds(i).neurons(inds);
        SummaryStruct_tmp.birds(1).birdname = SummaryStruct.birds(i).birdname;
        
        % === extract for just this expt
        [MOTIFSTATS] = lt_neural_v2_ANALY_ExtractMotif(SummaryStruct_tmp, ...
            collectWNhit, 0, LearnKeepOnlyBase);
        
        % === OUTPUT
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS = MOTIFSTATS;
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct = SummaryStruct_tmp;
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).exptname = exptname;
        MOTIFSTATS_Compiled.birds(i).birdname = birdname;

    end
end

