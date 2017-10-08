%% 9/28 -

%% =================== EXTRACT SOBER/WOHL RA DATA
clear all; close all;
SummaryStruct = lt_neural_RASamMel_SummaryStruct;


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
TimeWindowDur = 0.02;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-0.01:TimeWindowSlide:0.09; ...
    0.01:TimeWindowSlide:0.11]'; % N x 2 (pre and post onset)
% ListOfTimeWindows = [0.095:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
%     0.095+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.004];
savenotes = 'RAallbirds20ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


%%  AFP 20ms

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
TimeWindowDur = 0.02;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.004];
savenotes = 'LMANXallbirds20ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 1; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


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
