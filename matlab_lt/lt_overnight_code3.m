close all;

%% 1
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 2;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly);


% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);

% -- save
lt_save_figs_to_folder('X_includingLearning',1);


%% 2
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'X', 'LMAN'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly);


% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);

% -- save
lt_save_figs_to_folder('X_LMAN_includingLearning',1);



%% 2
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 2;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly);


% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);

% -- save
lt_save_figs_to_folder('LMAN_noLearning',1);


%% 2
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly);


% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);

% -- save
lt_save_figs_to_folder('LMAN_withLearning',0);


%% ========================== 11/18/15 - wh25 context autolabel
clear all; close all;
batch = 'batch.keep';
% config='/bluejay4/lucas/birds/wh25pk77/111515_CtxtDepPitch_away_durWN/CONFIG_away.evconfig2';
config='/bluejay4/lucas/birds/wh25pk77/config_112015.evconfig2';


syl.targ='b';
syl.pre='';
syl.post=''; 
NoteNum=0; 

ampThresh=4000;
min_dur=20;
min_int=4;

overwrite_notmat=1;

% ---- RUN
MetadataStruct=lt_metadata_collect;

experiment = 'CtxtDepPitch';
condition='';
notes='';
date_range={'06Dec2015', '14Dec2015'};
only_labeled_dirs=0;
ListOfDirs2=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs);

for i=1:length(ListOfDirs2);
   
    cd(ListOfDirs2{i});
   
   lt_make_batch(1);
    
[fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);

cd ..
end


