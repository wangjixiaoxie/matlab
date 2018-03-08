%% ###############################################################
%% ##################################################### POPULATION
%% ======================== EXTRACT SEGMENTS FOR POPULATIONS
% ========== KEY DIFFERENCES,
% 1) MIGHT HAVE TO SKIP SONGS TO MAKE SURE ALL
% NEURONS ARE REPRESENTED IN ALL SONGS
% 2)


close all; clear MOTIFSTATS_Compiled;
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 1;
saveOn = 1;
OrganizeByExpt =0;
collectFF=1;
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF);



%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
clear MOTIFSTATS_Compiled;


%% ================ PLOT [CORRELATION WITH FF]
close all;
MOTIFSTATS_pop = lt_neural_POP_FFcorr(MOTIFSTATS_pop, SummaryStruct);

close all;
lt_neural_POP_FFcorrPlot

%% ================ 

lt_neural_POP_PlotRast






%% #############################################################
%% ######################## POPULATION - TAKE ENTIRE MOTIF

% ================ extract full motifs and time warp here
lt_neural_v2_ExtractFullMotifs;


%% ======================== EXTRACT POPULATION
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
clear MOTIFSTATS_Compiled;


%% ================ PLOT [CORRELATION WITH FF]
close all;
xcov_dattotake = [];
xcovwindmax = 0.2;
binsize_spk = 0.005;
MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk);


%% =============== SUMMARY PLOT OF ALL CROSS-CORRELATIONS
close all;
OUTSTRUCT = lt_neural_POP_PlotSummary(MOTIFSTATS_pop, SummaryStruct);

lt_neural_POP_PlotSummary2(MOTIFSTATS_pop, SummaryStruct, OUTSTRUCT);


