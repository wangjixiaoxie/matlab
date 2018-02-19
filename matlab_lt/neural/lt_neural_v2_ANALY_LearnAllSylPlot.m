function lt_neural_v2_ANALY_LearnAllSylPlot(MOTIFSTATS_Compiled, ...
    MeanSubtract, BirdsToPlot)

% BirdsToPlot = {'pu69wh78'}, leave emxpty if plot all.
%% lt 11/5/17 - plots learning for all training expts (separates into targ, nontarg ctxt, nontarg syl)

NumBirds = length(MOTIFSTATS_Compiled.birds);
for i=1:NumBirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    
    if ~isempty(BirdsToPlot)
       if ~any(strcmp(BirdsToPlot, birdname))
           continue
       end
    end
%     ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
    ListOfExpts = unique({MOTIFSTATS_Compiled.birds(i).exptnum.exptname});
    for ll = 1:length(ListOfExpts)
        
        MOTIFSTATS = MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS;
        SummaryStruct_onebird = MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct;
        
%         
%         birdname = SummaryStruct.birds(i).birdname;
        exptname = ListOfExpts{ll};
        
        if ~isfield(MOTIFSTATS.params, 'TargSyls')
            disp('PROBLEM - TargSyls not in Params');
        end
        
        if isempty(MOTIFSTATS.params.TargSyls)
            sdafasdf
        end
        
        % =============== PLOT FOR THIS BIRD/EXPT
        lt_neural_v2_ANALY_LearnAllSylPlot2(SummaryStruct_onebird, MOTIFSTATS, ...
            MeanSubtract);
        
    end
end
