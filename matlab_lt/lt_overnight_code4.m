clear CLASSES

% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
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
ListOfFrBinSizes = [0.005 0.01];
savenotes = 'allXLman40ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.08;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005 0.01];
savenotes = 'allXLman80ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


% =============== 1
% ITERATIONS
TimeWindowDur = 0.03;
TimeWindowSlide = 0.005;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005];
savenotes = 'allXLman30ms';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);
