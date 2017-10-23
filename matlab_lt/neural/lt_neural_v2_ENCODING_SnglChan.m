function lt_neural_v2_ENCODING_SnglChan(SummaryStruct)
%% for each neuron calculate FF correlation with each syl (separated by context)


%% 
collectWNhit=0;
LearnKeepOnlyBase = 1; 
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, collectWNhit, LearnKeepOnlyBase);


%%