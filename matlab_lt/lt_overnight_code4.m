% analyfname = 'xaa_Algn2Ons1_27Oct2017_1114_XLMAN25ms';
% Niter = 1000;
% TimeWindows = [-0.035 -0.035]; % [-0.05 -0.05] means window from 50ms pre onset to 50ms pre offset (each row is separate analysis)
% lt_neural_v2_CTXT_BRANCH_DatVsShuff(analyfname, Niter, TimeWindows);
% 
% 
% 
% analyfname = 'xaa_Algn2Ons1_27Oct2017_1156_RA25ms';
% Niter = 1000;
% TimeWindows = [-0.020 -0.020]; % [-0.05 -0.05] means window from 50ms pre onset to 50ms pre offset (each row is separate analysis)
% lt_neural_v2_CTXT_BRANCH_DatVsShuff(analyfname, Niter, TimeWindows);


%%
% 
% clear all; close all;
% BirdsToKeep = {'bu77wh13'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'LMAN', 'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 0;
% BatchesDesired = {};
% ChannelsDesired = [];
% % BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% % BrainArea = {'X'};
% % ExptToKeep = {};
% % RecordingDepth = [];
% % LearningOnly = 1;
% % BatchesDesired = {};
% % ChannelsDesired = [];
% [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
%     BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);
% 
% % --- load all neurons
% if (0)
%     [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
% end
% 
% % ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM
% SummaryStruct =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);
% 
% 
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'xaa'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
% 
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
% prms.alignWhichSyl = 2; % which syl (in order) to align to
% prms.alignOnset = 1; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.15;
% prms.motifpostdur = 0.15;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
% 
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFIER (V2) - picks a branch,
% % does all time points, goes to next branch.
% TimeWindowDur = 0.04;
% TimeWindowSlide = 0.015;
% FRbinsize = 0.005;
% savenotes = 'testLMAN1birdsDiffGLMnet';
% 
% prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassSlide.GetPosControl =1;
% 
% CVmethod = 'Kfold';
% plotstat = 'F1';
% 
% saveON =1;
% 
% ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
%     TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
%     saveON);

%% xaaa

clear all; close all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN', 'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end

% ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM
SummaryStruct =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.15;
prms.motifpostdur = 0.35;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFIER (V2) - picks a branch,
% does all time points, goes to next branch.
TimeWindowDur = 0.025;
TimeWindowSlide = 0.005;
FRbinsize = 0.005;
savenotes = 'XLMAN25msLTW';

prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassSlide.GetPosControl =1;

CVmethod = 'Kfold';
plotstat = 'F1';

saveON =1;
LinTimeWarp = 1;
regionstowarp = [3 4 5 6];
ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
    TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
    saveON, LinTimeWarp, regionstowarp);



%% xaaa

clear all; close all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN', 'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end

% ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM
SummaryStruct =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 4; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.45;
prms.motifpostdur = 0.1;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFIER (V2) - picks a branch,
% does all time points, goes to next branch.
TimeWindowDur = 0.025;
TimeWindowSlide = 0.005;
FRbinsize = 0.005;
savenotes = 'XLMAN25msLTW';

prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassSlide.GetPosControl =1;

CVmethod = 'Kfold';
plotstat = 'F1';

saveON =1;
LinTimeWarp = 1;
regionstowarp = [3 4 5 6];
ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
    TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
    saveON, LinTimeWarp, regionstowarp);


%% xaaa

clear all; close all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN', 'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end

% ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM
SummaryStruct =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 4; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.45;
prms.motifpostdur = 0.1;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFIER (V2) - picks a branch,
% does all time points, goes to next branch.
TimeWindowDur = 0.025;
TimeWindowSlide = 0.005;
FRbinsize = 0.005;
savenotes = 'XLMAN25msLTW';

prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassSlide.GetPosControl =1;

CVmethod = 'Kfold';
plotstat = 'F1';

saveON =1;
LinTimeWarp = 1;
regionstowarp = [3 4 5 6 7 8];
ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
    TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
    saveON, LinTimeWarp, regionstowarp);


%% =================== EXTRACT SOBER/WOHL RA DATA
% clear all; close all;
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'RA'};
% % ExptToKeep = {'RAlearn1', 'RALMANlearn1', 'LMANsearch'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 0;
% BatchesDesired = {};
% ChannelsDesired = [];
% % BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% % BrainArea = {'X'};
% % ExptToKeep = {};
% % RecordingDepth = [];
% % LearningOnly = 1;
% % BatchesDesired = {};
% % ChannelsDesired = [];
% [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
%     BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);
% 
% 
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'xaa'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
% 
% 
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
% prms.alignWhichSyl = 2; % which syl (in order) to align to
% prms.alignOnset = 1; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.15;
% prms.motifpostdur = 0.15;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
% 
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFIER (V2) - picks a branch,
% % does all time points, goes to next branch.
% TimeWindowDur = 0.025;
% TimeWindowSlide = 0.005;
% FRbinsize = 0.005;
% savenotes = 'RA25msLTW';
% 
% prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassSlide.GetPosControl =1;
% 
% CVmethod = 'Kfold';
% plotstat = 'F1';
% 
% saveON =1;
% LinTimeWarp = 1;
% regionstowarp = [3 4];
% ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
%     TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
%     saveON, LinTimeWarp, regionstowarp);
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% continue where left off
% clear all; close all;
% load('SUMMARYv2_xaa_Algn2Ons1_01Dec2017_1653_RA25msLTW.mat');
% load('PARAMSv2_xaa_Algn2Ons1_01Dec2017_1653_RA25msLTW.mat');
% load('ALLBRANCHv2_xaa_Algn2Ons1_01Dec2017_1653_RA25msLTW.mat');
% load('CLASSESv2_xaa_Algn2Ons1_01Dec2017_1653_RA25msLTW.mat');
% tstampsave = '01Dec2017_1653';
% 
% TimeWindowDur = 0.025;
% TimeWindowSlide = 0.005;
% FRbinsize = 0.005;
% savenotes = 'RA25msLTW';
% 
% prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassSlide.GetPosControl =1;
% 
% CVmethod = 'Kfold';
% plotstat = 'F1';
% 
% saveON =1;
% LinTimeWarp = 1;
% regionstowarp = [3 4];
% ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
%     TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
%     saveON, LinTimeWarp, regionstowarp, ALLBRANCH, tstampsave);
% 
% 
% 

%% #############################################################
%% 9/28 -

%% =================== EXTRACT SOBER/WOHL RA DATA
% clear all; close all;
% SummaryStruct = lt_neural_RASamMel_SummaryStruct;
% 
% 
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'xaa'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
% 
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
% prms.alignWhichSyl = 2; % which syl (in order) to align to
% prms.alignOnset = 1; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.15;
% prms.motifpostdur = 0.2;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
% 
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
% 
% 
% % =============== 1
% % ITERATIONS
% TimeWindowDur = 0.02;
% TimeWindowSlide = 0.005;
% ListOfTimeWindows = [-0.01:TimeWindowSlide:0.09; ...
%     0.01:TimeWindowSlide:0.11]'; % N x 2 (pre and post onset)
% % ListOfTimeWindows = [0.095:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
% %     0.095+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
% ListOfFrBinSizes = [0.004];
% savenotes = 'RAallbirds20ms';
% 
% prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
% prms.ClassGeneral.GetPosControl =1;
% 
% [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
%     ListOfFrBinSizes, savenotes);
% 

%%  AFP 30ms

clear all; close all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN', 'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end


% ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM
SummaryStruct =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.15;
prms.motifpostdur = 0.2;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.03;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'LMANXallbirds30ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


%%
close all;
clear all;

strtype = 'xaa';
algnsyl = 2;
algnonset = 1;
suffix = 'LMANXallbirds20ms'; % leave blank to get any
CLASSEScompiled = lt_neural_v2_CTXT_PlotGeneral_M(strtype, algnsyl, ...
    algnonset, suffix);


close all; 
strtype = 'xaa';
plotstat = 'F1';
suffix = 'LMANXallbirds20ms'; % leave blank to get all
ALLBRANCH = lt_neural_v2_CTXT_PlotAll(strtype, plotstat, suffix);


%% ########################### below, night of 9/27
%% RA data - [xaaa, syl3]

% =================== EXTRACT SOBER/WOHL RA DATA
clear all; close all;
SummaryStruct = lt_neural_RASamMel_SummaryStruct;


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 3; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.3;
prms.motifpostdur = 0.2;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'RAallbirds';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);



% ####################################################### RA data - [xaaa, syl2]

% =================== EXTRACT SOBER/WOHL RA DATA
clear all; close all;
SummaryStruct = lt_neural_RASamMel_SummaryStruct;


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.18;
prms.motifpostdur = 0.2;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfTimeWindows = [-0.04:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -0.04+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'RAallbirds';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


%% LMAN syl 2

% EXTRACT 
clear all; close all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end


% ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM
SummaryStruct =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);



% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'xaa'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
% 
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
% prms.alignWhichSyl = 2; % which syl (in order) to align to
% prms.alignOnset = 1; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.17;
% prms.motifpostdur = 0.2;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
% 
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
% 
% 
% % =============== 1
% % ITERATIONS
% TimeWindowDur = 0.04;
% TimeWindowSlide = 0.01;
% ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
%     -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
% ListOfFrBinSizes = [0.005];
% savenotes = 'LMANallbirds';
% 
% prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
% prms.ClassGeneral.GetPosControl =1;
% 
% [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
%     ListOfFrBinSizes, savenotes);




% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.14;
prms.motifpostdur = 0.28;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'LMANallbirds';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);



%% RA data - [xaaa, syl3]

% =================== EXTRACT SOBER/WOHL RA DATA
clear all; close all;
SummaryStruct = lt_neural_RASamMel_SummaryStruct;


% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 3; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.3;
prms.motifpostdur = 0.2;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
CVmethod = 'Kfold';
CVkfoldnum = 8;
TimeWindowDur = 0.04;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'RAallbirds';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);






%% % %%  get all compiled
% 
% listofresults = dir('Results_*');
% 
% for i=1:length(listofresults)
%     uscores = strfind(listofresults(i).name, '_');
%     
%     strtype = listofresults(i).name(uscores(1)+1:uscores(2)-1);
%     algnsyl = listofresults(i).name(uscores(2)+8);
%     algnonset = listofresults(i).name(uscores(3)-1);
%     assert(uscores(3)-uscores(2) == 15, 'problem, mult digits')
%     
%     CLASSEScompiled = lt_neural_v2_CTXT_PlotGeneral_M(strtype, algnsyl, algnonset);
% end
% 
% 
% %%  then get all branches
%     cd('/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M')
% 
% listofresults = dir('CLASSEScompiled_*');
% 
% StrtypeAll = {};
% 
% for i=1:length(listofresults)
%     
%     uscores = strfind(listofresults(i).name, '_');
%     
%     strtype = listofresults(i).name(uscores(1)+1:uscores(2)-1);
%     StrtypeAll = [StrtypeAll strtype];
% end
% 
% StrtypeAll = unique(StrtypeAll);
% 
% for i=1:length(StrtypeAll)
%    
%     strtype = StrtypeAll{i};
%     close all;
%     plotstat = 'F1';
%     ALLBRANCH = lt_neural_v2_CTXT_PlotAll(strtype, plotstat);
%     
%     % ================ EXTRACT NEURAL FR TO BRANCHES
%     saveOn = 1;
%     ALLBRANCH = lt_neural_v2_CTXT_BranchGetFR(ALLBRANCH, saveOn);
%     
%     clear ALLBRANCH
%         cd('/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M')
% 
% end
% 

%% aax (2,1)
% clear CLASSES
%
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'aax'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
%
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
% prms.alignWhichSyl = 2; % which syl (in order) to align to
% prms.alignOnset = 1; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.2;
% prms.motifpostdur = 0.3;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
%
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
%
%
% % =============== 1
% % ITERATIONS
% TimeWindowDur = 0.04;
% TimeWindowSlide = 0.01;
% ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
%     -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
% ListOfFrBinSizes = [0.005];
% savenotes = 'allXLman40ms';
%
% prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
% prms.ClassGeneral.GetPosControl =1;
%
% [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
%     ListOfFrBinSizes, savenotes);
%
%
%% aax (1, 0)
% clear CLASSES
%
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'aax'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
%
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
% prms.alignWhichSyl = 1; % which syl (in order) to align to
% prms.alignOnset = 0; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.1;
% prms.motifpostdur = 0.4;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
%
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
%
%
% % =============== 1
% % ITERATIONS
% TimeWindowDur = 0.04;
% TimeWindowSlide = 0.01;
% ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
%     -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
% ListOfFrBinSizes = [0.005];
% savenotes = 'allXLman40ms';
%
% prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
% prms.ClassGeneral.GetPosControl =1;
%
% [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
%     ListOfFrBinSizes, savenotes);
%
%
% %% aax (1, 0)
% clear CLASSES
%
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'xaa'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
%
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
% prms.alignWhichSyl = 3; % which syl (in order) to align to
% prms.alignOnset = 1; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.35;
% prms.motifpostdur = 0.15;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
%
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
%
%
% % =============== 1
% % ITERATIONS
% TimeWindowDur = 0.04;
% TimeWindowSlide = 0.01;
% ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
%     -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
% ListOfFrBinSizes = [0.005];
% savenotes = 'allXLman40ms';
%
% prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
% prms.ClassGeneral.GetPosControl =1;
%
% [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
%     ListOfFrBinSizes, savenotes);
%
% %% aaax (1,0)
% clear CLASSES
% 
% % &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
% strtype = 'aaax'; % a is fixed, x variable, across contexts
% [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);
% 
% % &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
% prms.alignWhichSyl = 1; % which syl (in order) to align to
% prms.alignOnset = 0; % if 1, then onset, if 0, then offset
% prms.motifpredur = 0.15;
% prms.motifpostdur = 0.4;
% prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
% CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);
% 
% % &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
% CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);
% 
% 
% % =============== 1
% % ITERATIONS
% TimeWindowDur = 0.04;
% TimeWindowSlide = 0.01;
% ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
%     -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
% ListOfFrBinSizes = [0.005];
% savenotes = 'allXLman40ms';
% 
% prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
% prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
% prms.ClassGeneral.GetPosControl =1;
% 
% [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
%     ListOfFrBinSizes, savenotes);
% 

%% aaax (2,0)
clear CLASSES

% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'aaax'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.15;
prms.motifpostdur = 0.4;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'allXLman40ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);

%% aaax (3,0)
clear CLASSES

% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'aaax'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
prms.alignWhichSyl = 3; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.3;
prms.motifpostdur = 0.25;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'allXLman40ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);

%% xaaa (3,1)
clear CLASSES

% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
prms.alignWhichSyl = 3; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.3;
prms.motifpostdur = 0.25;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'allXLman40ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


%% xaaa (3,1)
clear CLASSES

% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 0; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.2;
prms.motifpostdur = 0.3;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'allXLman40ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);
