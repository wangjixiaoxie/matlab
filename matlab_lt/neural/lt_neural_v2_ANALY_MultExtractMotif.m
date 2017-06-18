function MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct)


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
        [MOTIFSTATS] = lt_neural_v2_ANALY_ExtractMotif(SummaryStruct_tmp);
        
        % === OUTPUT
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS = MOTIFSTATS;
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct = SummaryStruct_tmp;
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).exptname = exptname;
        MOTIFSTATS_Compiled.birds(i).birdname = birdname;

    end
end
