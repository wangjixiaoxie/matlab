function lt_neural_v2_ANALY_LearnAllSylPlot(SummaryStruct, MOTIFSTATS_Compiled, ...
    MeanSubtract)

%% lt 11/5/17 - plots learning for all training expts (separates into targ, nontarg ctxt, nontarg syl)

close all;
NumBirds = length(SummaryStruct.birds);
for i=1:NumBirds
    ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
    
    for ll = 1:length(ListOfExpts)
        
        MOTIFSTATS = MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS;
        SummaryStruct_onebird = MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct;
        
        birdname = SummaryStruct.birds(i).birdname;
        exptname = ListOfExpts{ll};
        
        if ~isfield(MOTIFSTATS.params, 'TargSyls')
            disp('PROBLEM - TargSyls not in Params');
        end
        if isempty(MOTIFSTATS.params.TargSyls);
            sdafasdf
        end
        
        % =============== PLOT FOR THIS BIRD/EXPT
        lt_neural_v2_ANALY_LearnAllSylPlot2(SummaryStruct_onebird, MOTIFSTATS, ...
            MeanSubtract);
        

        
        
    end
end
