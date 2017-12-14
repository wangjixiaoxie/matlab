clear all; close all;
% =================== HERE LOAD THE PARTIALLY COMPLETE STRUCTURES
load('SUMMARYv2_xaa_Algn2Ons1_30Nov2017_1645_RA25msLTW.mat');
load('PARAMSv2_xaa_Algn2Ons1_30Nov2017_1645_RA25msLTW.mat');
load('ALLBRANCHv2_xaa_Algn2Ons1_30Nov2017_1645_RA25msLTW.mat');
load('CLASSESv2_xaa_Algn2Ons1_30Nov2017_1645_RA25msLTW.mat');
tstampsave = '30Nov2017_1645';

% ================== THE REST PARAMS ENTER AS BEFORE
TimeWindowDur = 0.025;
TimeWindowSlide = 0.005;
FRbinsize = 0.005;
savenotes = 'RA25msLTW';

prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassSlide.GetPosControl =1;

CVmethod = 'Kfold';
plotstat = 'F1';

saveON =1;
LinTimeWarp = 1;
regionstowarp = [3 4];
ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
    TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
    saveON, LinTimeWarp, regionstowarp, ALLBRANCH, tstampsave);

