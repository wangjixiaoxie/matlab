%% ACROSS NEURONS, FOR ANALYSES
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% EXTRACT 

% BirdsToKeep = {'bk7', 'wh6pk36', 'br92br54', 'or74bk35'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'LMAN', 'X'};
% ExptToKeep = {'LearnLMAN1', 'LMANlearn2', 'LMANneural2'};
% RecordingDepth = [1860 2770 2840 1800 2490];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'X', 'LMAN'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 1;
BatchesDesired = {};
ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end

%% ==== refinalize specific neurons, with slight changes
if (0)
birdnum=1;
i=6;

cd(SummaryStruct.birds(birdnum).neurons(i).dirname);
SummaryStruct.birds(birdnum).neurons(i)

clustnum = SummaryStruct.birds(birdnum).neurons(i).clustnum;
depth = SummaryStruct.birds(birdnum).neurons(i).electrode_depth;
Notes = SummaryStruct.birds(birdnum).neurons(i).Notes;
LEARN_WNonDatestr = SummaryStruct.birds(birdnum).neurons(i).LEARN_WNonDatestr;
LEARN_WNotherImportantDates = SummaryStruct.birds(birdnum).neurons(i).LEARN_WNotherImportantDates;

% CHANGES
Notes{2} = 'Location_notLMAN';

% RUN
lt_neural_v2_Finalize(clustnum, depth, Notes, {LEARN_WNonDatestr, LEARN_WNotherImportantDates})
end

%% EXTRACT FF AND SAVE

FFparamsAll.bird(1).birdname = 'bk7';
FFparamsAll.bird(1).FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
            'v', [2450 4300]};
FFparamsAll.bird(1).FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
            'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
FFparamsAll.bird(1).FFparams.cell_of_FFtimebins_DurLearn={'h', [0.034 0.038], 'b', [0.053 0.07], ...
            'v', [0.052 0.07]}; % WN on g H

       
FFparamsAll.bird(2).birdname = 'bu77wh13';
FFparamsAll.bird(2).FFparams.cell_of_freqwinds={'b', [2700 3900], 'h', [2600 3900], 'a', [1300 2600]};
FFparamsAll.bird(2).FFparams.cell_of_FFtimebins={'b', [0.031 0.04], 'h', [0.038 0.052], 'a', [0.056 0.081]}; % in sec, relative to onset (i.e. see vector T)
FFparamsAll.bird(2).FFparams.cell_of_FFtimebins_DurLearn={'b', [0.031 0.04], 'h', [0.038 0.052], 'a', [0.056 0.081]}; % in sec, relative to onset (i.e. see vector T)


FFparamsAll.bird(3).birdname = 'wh6pk36';
FFparamsAll.bird(3).FFparams.cell_of_freqwinds={'c', [2100 3100], 'h', [2800 4000], 'b', [2700 3800], ...
    'a', [1300 2200], 's', [4000 5100], 'd', [900 2000],  'n', [3300 4300], 'v', [2600 4000]};
FFparamsAll.bird(3).FFparams.cell_of_FFtimebins={'c', [0.035 0.05], 'h', [0.023 0.033], 'b', [0.034 0.038], ...
    'a', [0.06 0.08], 's', [0.02 0.023], 'd', [0.02 0.027],  'n', [0.029 0.05], 'v', [0.028 0.036]};
FFparamsAll.bird(3).FFparams.cell_of_FFtimebins_DurLearn={'c', [0.035 0.05], 'h', [0.023 0.033], 'b', [0.034 0.038], ...
    'a', [0.06 0.08], 's', [0.02 0.023], 'd', [0.02 0.027],  'n', [0.029 0.05], 'v', [0.028 0.036]};


% TEMPORARY
FFparamsAll.bird(4).birdname = 'br92br54';
FFparamsAll.bird(4).FFparams.cell_of_freqwinds={'a', [750 1400], 'c', [1200 1800], ...
    'h', [2350 3900], 'd', [1300 3400], 'k', [800 1800]};
FFparamsAll.bird(4).FFparams.cell_of_FFtimebins={'a', [0.057 0.08], 'c', [0.044 0.055], ...
    'h', [0.040 0.049], 'd', [0.032 0.052], 'k', [0.05 0.055]};
FFparamsAll.bird(4).FFparams.cell_of_FFtimebins_DurLearn={'a', [0.057 0.08], 'c', [0.044 0.055], ...
    'h', [0.040 0.049], 'd', [0.032 0.052], 'k', [0.05 0.055]};


FFparamsAll.bird(5).birdname = 'or74bk35';
FFparamsAll.bird(5).FFparams.cell_of_freqwinds={'b', [2750 3900]};
FFparamsAll.bird(5).FFparams.cell_of_FFtimebins={'b', [0.033 0.041]};
FFparamsAll.bird(5).FFparams.cell_of_FFtimebins_DurLearn={'b', [0.033 0.041]};


overWrite = 0; % note, will overwrite rgardless if detects chagnes (NOTE: always overwrites if detects changes)
plotSpec = 0; % to plot raw spec overlayed with PC and windows.
plotOnSong = 25; % will only start plotting spec once hit this song num.
plotSyl = ''; % to focus on just one syl. NOT DONE YET
lt_neural_v2_EXTRACT_FF_tmp(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl);


%% ======== SAVE META INFO FOR LEARNING EXPT HERE
% edit this by hand

% 1) --- LOAD
LearningMetaDat = lt_neural_v2_LoadLearnMetadat;


% 2) --- EDIT
LearningMetaDat; % OPEN AND EDIT BY HAND. 
% Note: each expt/targ syl has one column:
% row 3 and larger gives time of switch and nature of switch (4 types of
% switches possible) (escape dir)
% Al = 100% WN
% Of = off
% Up = escape up
% Dn = escape dn


% 3) --- SAVE
cd('/bluejay5/lucas/analyses/neural/');
save('LearningMetaDat', 'LearningMetaDat');

%% ==== LIST OF MOTIFS FOR EACH BIRD/EXPERIMENT
if (0) % RUNS AUTOMATICALLY WHEN EXTRACT
SummaryStruct = lt_neural_v2_PostInfo(SummaryStruct);
end
%% ===== EXTRACT WN HITS

% IN PROGRESS ! see inside
lt_neural_v2_EXTRACT_WNhit(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl);



%% ======== VOCLATIONS --> predict neural

close all;
plotRaw = 0; % individaul neuronf igs (note, hand entering neuron ind currently, in code)
% Binparams.Pretime = 0.2; % to start getting binned data (rel to onset)
% Binparams.Posttime = 0.2; % to stop getting dat (rel to onset);
% Binparams.Binsize = 0.025; % for getting spike counts
Binparams.Pretime = 0.7; % to start getting binned data (rel to onset)
Binparams.Posttime = 0.7; % to stop getting dat (rel to onset);
Binparams.Binsize = 0.025; % for getting spike counts

ConvOrDiv = 'conv';
saveOn = 0;

VOCALSTRUCTall = lt_neural_v2_ANALY_VocModel(SummaryStruct, Binparams, plotRaw, ConvOrDiv, saveOn);


% ===== 1) FILTER VOCAL STRUCT


% ===== 2) PLOT
close all;
birdnum = 1;
neuronnum=1;
timebin = 35;
plotSmoothed = 0; 
lt_neural_v2_ANALY_VocModel_plot(VOCALSTRUCTall, birdnum, neuronnum, timebin, plotSmoothed)


% ===== 3) ANOVA (done for each time bin)
% currerntly good. apply filter preceding to extract,e.g.. only one branch
% point, and so on. 
lt_neural_v2_ANALY_VocModel_anova(VOCALSTRUCTall)


% ====== 4) 
lt_neural_v2_ANALY_VocModel_glm(VOCALSTRUCTall)



%% ===== EFFECT OF MOTIF POSITION (WITHIN SONG BOUT) ON SYL FEATURES + NEURAL [JOANNE]
% INSERT INTO SUMMARY STRUCT, DON'T NEED TO SAVE AS IS VERY QUICK TO RUN
% THIS EVERY TIME
close all;

BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly);


% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
close all;
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++ LEARNING

%% ============== FOR EACH LEARNING EXPERIMENT, PLOT ALL NEURONS AND ALL MOTIFS
% ONE BIRD/LEARNING EXPT AT A TIME [BUT CAN COMBINE NEURONS]

% === 1) CHOOSE ONE LEARNING EXPT
% ===== LMANlearn2
BirdsToKeep = {'br92br54', 'or74bk35'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {};
ExptToKeep = {'LMANlearn4', 'LMANneural2'};
RecordingDepth = [1975 2490];
LearningOnly = 1; % then only if has WN on/switch time date
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly);

% ==== 1) EXTRACT DATA FOR EACH NEURON, EACH MOTIF
[MOTIFSTATS, SummaryStruct_filt] = lt_neural_v2_ANALY_LearningExtractMotif(SummaryStruct);

close all;
% motifs for this bird
% learning expt id
plottype = 'byneuron'; % 
% 'byneuron' - each neuron one fig will all motifs [DEFAULT]
% 'dotprod' - for each bin of trials get dot prod from IN PROGRESS
% 'bysyl' - each plot one syl, all neurons.
DivisorBaseSongs = 1;
lt_neural_v2_ANALY_Learning(SummaryStruct_filt, MOTIFSTATS, plottype, DivisorBaseSongs);

% ===== FOR ALL NEURONS, PLOT CHANGE AT TARG VS. OTHERS. ALSO, AVERAGE
% DEVIATIONS ACROSS ALL EXPERIMENTS.


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% =================== LEARNING, SUMMARY STATISTICS ACROSS ALL EXPERIMENTS
% CAN DO MULTIPLE BIRDS

NumBirds = length(SummaryStruct.birds);
MOTIFSTATS_Compiled = struct; % heirarchy: birds --> expt --> neurons

for i=1:NumBirds

    ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
    
    for ll=1:length(ListOfExpts)
        exptname = ListOfExpts{ll};
        
        inds = strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname);
        
        SummaryStruct_tmp = struct;
        SummaryStruct_tmp.birds(1).neurons = SummaryStruct.birds(i).neurons(inds);
        SummaryStruct_tmp.birds(1).birdname = SummaryStruct.birds(i).birdname;
        
        % === extract for just this expt
        [MOTIFSTATS] = lt_neural_v2_ANALY_LearningExtractMotif(SummaryStruct_tmp);
        
        % === OUTPUT
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS = MOTIFSTATS;
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct = SummaryStruct_tmp;
        MOTIFSTATS_Compiled.birds(i).exptnum(ll).exptname = exptname;

    end
end


%% =========== PICK OUT SAME TYPE/DIFF

numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:numbirds
    
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
   
    for ii=1:numexpts
       
        TargSyls = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.TargSyls;
        MotifsActual = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.MotifsActual;
        
        [SameTypeSyls, DiffTypeSyls, motif_regexpr_str] = ...
        lt_neural_v2_extractSameType(MotifsActual, TargSyls);

    % --- OUTPUT
MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.SameTypeSyls = SameTypeSyls;
MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.DiffTypeSyls = DiffTypeSyls;


    end
end


%% ==== FOR EACH LEARNING EXPERIMENT, PLOT TIMELINE OF NEURONS AND LEARNING
for i=1:NumBirds
   ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
   
   for ll = 1:length(ListOfExpts)
        
       MOTIFSTATS = MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS;
       SummaryStruct_tmp = MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct;
       
      % === PLOT
      lt_neural_v2_ANALY_LearningPlot1(SummaryStruct_tmp, MOTIFSTATS);
       
   end
end
    
%% ==== FOR EACH LEARNING EXPERIMENT, PLOT SUMMARY STATISTICS
close all;

PlotTimeCourses = 0; % plots raw timecourse for each expt/neuron
PlotNeuralFFDeviationCorrs = 0; % zscore of neural similarity and FF, correalted fore aach syl
convertToZscore = 0;

LearningMetaDat = lt_neural_v2_LoadLearnMetadat;
lt_neural_v2_ANALY_LearningStats(MOTIFSTATS_Compiled, LearningMetaDat, PlotTimeCourses, ...
    PlotNeuralFFDeviationCorrs, convertToZscore);







