%% ###########################################################################
%% ############################################## DATA PREPROCESSING
% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
clear CLASSES

strtype = 'xaaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.15;
prms.motifpostdur = 0.1;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% &&&&&&&&&&&&& 2) PLOT MEAN FR ACROSS CONTEXTS FOR EACH BRANCH 
close all;
plotPosControl = 0; % will do if exists.
LMANorX = 1; % 0 for all; 1 for LMAN; 2 for X
closeAfterEachBird = 1; % closes figss
lt_neural_v2_CTXT_FRanyclass(CLASSES, SummaryStruct, prms, plotPosControl, ...
    LMANorX, closeAfterEachBird);


%% ###########################################################################
%% ##################   [OLD VERSION]
%% ================================ CLASSIFICATION (SINGLE TIME WINDOW)

% ========================= CLASSIFIER( SINGLE ITERATION)
prms.ClassGeneral.frtimewindow =[-0.075 0.025]; % on and off, relative to syl onset
prms.ClassGeneral.frbinsize = 0.01; % in s.
prms.ClassGeneral.Nmin = 7; % in s.

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case

prms.ClassGeneral.GetPosControl =1;
CLASSES = lt_neural_v2_CTXT_ClassGeneral(CLASSES, SummaryStruct, prms);


% ===================== PLOTS SINGLE ITERATION
close all;
lt_neural_v2_CTXT_PlotGeneral(CLASSES, SummaryStruct, prms);



%% ================================ CLASSIFICATION (SLIDING TIME WINDOW)


% ============================ CLASSIFIER (SLIDING WINDOW MULT ITERATIONS)
TimeWindowDur = 0.04;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005 0.01 0.02];
savenotes = 'allXLman';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case
prms.ClassGeneral.GetPosControl =1;

[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);


%% ========================== COMPILING AND PLOTTING RESULTS
% ======================================== PLOT RESULTS FOR SINGLE ANALYSIS
% 1) COLLECT DATA
strtype = 'xaa';
algnsyl = 2;
algnonset = 1;
suffix = 'RAallbirds20ms'; % leave blank to get any
CLASSEScompiled = lt_neural_v2_CTXT_PlotGeneral_M(strtype, algnsyl, ...
    algnonset, suffix);

% 2) PLOT
close all; 
plotstat = 'F1';
lt_neural_v2_CTXT_PlotGeneral_M2(CLASSEScompiled, plotstat);



% ===================================== COMBINE ALL PLOTS FOR A GIVEN MOTIF 
% NOTE!!: NEED TO FIRST RUN lt_neural_v2_CTXT_PlotGeneral_M to get compiled
% stats
% combine across 1) syl aligned to 2) aligned to onset or offset, 3) fr
% window size 4) fr bin size
close all; 
strtype = 'xaa';
plotstat = 'F1';
suffix = 'RAallbirds20ms'; % leave blank to get all
ALLBRANCH = lt_neural_v2_CTXT_PlotAll(strtype, plotstat, suffix);


%% ###########################################################################
%% ##################   [NEW VERSION]
%% ===========================  CLASSIFIER (V2) - picks a branch,
% does all time points, goes to next branch.
TimeWindowDur = 0.025;
TimeWindowSlide = 0.005;
FRbinsize = 0.005;
savenotes = 'pu69wh78RALMAN';

prms.ClassSlide.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassSlide.GetPosControl =1;

CVmethod = 'Kfold';
plotstat = 'F1';

saveON =1;
LinTimeWarp = 1;
regionstowarp = [3 4];

ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
    TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
    saveON, LinTimeWarp, regionstowarp);

% NOTE: TO CONTINUE WHERE LEFT OFF BEFORE, RUN ABOVE FUNCTION WITH MODIFIED
% VARIABLES (SEE THIS SCRIPT FOR EXAMPLE)
lt_neural_v2_CTXT_ClassSliding_Continue;


% ------------ debugging: to systematically change names of classes...
lt_neural_v2_CTXT_Debug;




%% ############################# PLOTTING ALLBRANCh [SINGLE ANALYSIS]
% ======= EXTRACT GAP/SYL DURS
close all;
ALLBRANCH = lt_neural_v2_CTXT_BranchGaps(ALLBRANCH);

% ==== 1)  REMOVE ANY REDUNDANT NEURONS FROM ALLBRANCH
ALLBRANCH = lt_neural_v2_CTXT_BranchRemvOlap(ALLBRANCH);


% ==== 2)  PLOT EACH BRANCH/BIRD/NEURON
close all;
birdtoplot = 'pu69wh78'; % leave blank to plot all;
plotspec_num = 0; % how many spectrograms to plot for each class in each branch point? if 0 then none.
locationtoplot = {};
BranchToPlot = {'[a-z]bh'}; % type regexp strings
plotrasters = 0;
lt_neural_v2_CTXT_BranchEachPlot(ALLBRANCH, birdtoplot, plotspec_num, ...
    locationtoplot, BranchToPlot, plotrasters)


% ==== 3)  SUMMARIZE PLOT ACROSS BRANCHES 
close all; 
dattoplot = 'classperform';
% dattoplot = 'frmean';
% dattoplot = 'dprime';
LMANorX = 2; % 0, both; 1, LMAN; 2, X, 3(RA)
birdstoexclude = {};
% birdstoexclude = {'bk7', 'bu77wh13', 'or74bk35', 'wh6pk36', 'br92br54'};

% durThreshOmega.syl = 0.15; % omega2 (will only keep if lower) [leave empty to ignore]
% durThreshOmega.gappre= 0.5;
% durThreshOmega.gappost= 0.2;
durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
durThreshOmega.gappre= [];
durThreshOmega.gappost= [];

RemoveRepeats=0; % if 1, then removes any branch with a class with token preceded by same syl (e.g. a(a)bc or a(a)ab);
RemovePrecededByIntro=0; % if 1, then removes any token preceded by "i" (e.g. i(a)bc) [REQUIRES REMOVE REPEATS =1, since in progress]

lt_neural_v2_CTXT_PlotAllBranch(ALLBRANCH, LMANorX, dattoplot, birdstoexclude, ...
    durThreshOmega, RemoveRepeats, RemovePrecededByIntro)


% ################### [IMPORTANT] SORT BY BRANCH ID, BRAIN REGION, ETC
close all;
% BrainRegions = {'LMAN', 'X', 'RA'};
BrainRegions = {'LMAN', 'RA'};
% BrainRegions = {'LMAN','X'};
BirdToPlot = {};
useDprime=0;
lt_neural_v2_CTXT_BRANCH_PlotByBranchID(ALLBRANCH, BrainRegions, ...
    BirdToPlot, useDprime)


% ###################### PLOT EXAMPLES FOR EACH BIRD/BRANCH POINT




%% ############################# PLOTTING ALLBRANCH [MULTIPLE ANALYSIS]
% ################################## COMPARE TIMING OF TWO COMPILED BRANCHES
close all;
% branchfname1 = 'ALLBRANCH_xaa_03Oct2017_2149_RAallbirds20ms.mat';
% branchfname2 = 'ALLBRANCH_xaa_18Oct2017_2355_LMANXallbirds20ms.mat';
branchfname1 = 'ALLBRANCHv2_xaa_Algn2Ons1_30Nov2017_1911_XLMAN25msLTW.mat';
branchfname2 = 'ALLBRANCHv2_xaa_Algn2Ons1_01Dec2017_1653_RA25msLTW.mat';
% branchfname2 = 'ALLBRANCHv2_xaa_Algn2Ons1_29Nov2017_1838_pu69wh78RALMAN';
lt_neural_v2_CTXT_BranchCompareTwo(branchfname1, branchfname2);


% ====== TEMP, MOVE TO FUNCTION - SAVES COMPILED FOR ALL


%% ################################################### 
%% #################### [PREMOTOR WINDOW, DECODING] 
% ============= 1) IN PREMOTOR WINDOW, COMPARE DECODING VS. SHUFFLED.
close all;
analyfname = 'xaaa_Algn2Ons1_19Dec2017_1219_XLMAN25msLTW';
Niter = 1000;
TimeWindows = [-0.035 -0.035]; % [-0.05 -0.05] means window from 50ms pre onset to 50ms pre offset (each row is separate analysis)
% TimeWindows = [0 0 ]; % [-0.05 -0.05] means window from 50ms pre onset to 50ms pre offset (each row is separate analysis)
% TimeWindows = [-0.035 -0.035]; % LMAN
% TimeWindows = [-0.02 -0.02]; % RA
lt_neural_v2_CTXT_BRANCH_DatVsShuff(analyfname, Niter, TimeWindows);

% ------- to plot results from above (can do multiple)
close all;
% allanalyfnames = {...
%     'xaa_Algn2Ons1_30Nov2017_1911_XLMAN25msLTW', ...
%     };
allanalyfnames = {...
    'xaaa_Algn3Ons1_15Dec2017_0110_XLMAN25msLTW'};
allanalyfnames = {...
    'xaaa_Algn2Ons1_19Dec2017_1219_XLMAN25msLTW', ...
    'xaaa_Algn3Ons1_15Dec2017_0110_XLMAN25msLTW', ...
    'xaaa_Algn4Ons1_15Dec2017_1100_XLMAN25msLTW'};
DecodeStruct = lt_neural_v2_CTXT_BRANCH_DatVsShuffMULT(allanalyfnames);


%% ################ RUNNING HISTOGRAM DISTANCE AS DISTANCE METRIC

lt_neural_v2_CTXT_Distance();

numbirds = length(CLASSES.birds);
for i=1:numbirds
    
    numneurons = length(CLASSES.birds(i).neurons);
    
    for ii = 1:numneurons
       
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for iii=1:numbranches
           
            %%
                        if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT')
                continue
            end
            
            NN = [];
            for j=1:length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum)
                n = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(j).SegmentsExtract);
                NN = [NN n];
            end
            
            numrows = sum(NN >= prms.ClassSlide.Nmin);
            if numrows<2
                disp(['N too small!! - Skipped, ' birdname '-n' num2str(ii) '-br' num2str(iii)]);
                continue
            else
                disp(['N fine!! - continu, ' birdname '-n' num2str(ii) '-br' num2str(iii)]);
            end
            
            % *********************** confirm that pos control has same sampel size
            NN2 = [];
            for j=1:length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum)
                n = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(j).SegmentsExtract);
                NN2 = [NN2 n];
            end
            
            NN = NN(NN>=prms.ClassSlide.Nmin);
            NN2 = NN2(NN2>=prms.ClassSlide.Nmin);
            
            tmp = tabulate(NN)==tabulate(NN2);
            assert(all(tmp(:)==1), 'pos control not same sampoel size as dat');
            % ***********************
            
            ConfMatAll = nan(numrows, numrows, size(ListOfTimeWindows,1));
            ConfMatAll_NEG = nan(numrows, numrows, size(ListOfTimeWindows,1));
            ConfMatAll_POS = nan(numrows, numrows, size(ListOfTimeWindows,1));
            
            SEGEXTRACT_DAT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT;
            
            if LinTimeWarp ==1
                % ==== 1) combine all classes into one segextract
                segextractAll = SEGEXTRACT_DAT.classnum(1).SegmentsExtract;
                originalinds = [length(SEGEXTRACT_DAT.classnum(1).SegmentsExtract)];
                for j=2:length(SEGEXTRACT_DAT.classnum)
                    segextractAll = [segextractAll ...
                        SEGEXTRACT_DAT.classnum(j).SegmentsExtract];
                    
                    originalinds = [originalinds ...
                        length(SEGEXTRACT_DAT.classnum(j).SegmentsExtract)];
                end
                
                % ==== 2) linear time warp
                segextractAll = lt_neural_LinTimeWarpSegmented(segextractAll, regionstowarp);
                
                % ===== put back into SEGEXTRACT
                count=1;
                for j=1:length(SEGEXTRACT_DAT.classnum)
                    %                     disp(count:(count+originalinds(j)-1));
                    %                     disp('--');
                    SEGEXTRACT_DAT.classnum(j).SegmentsExtract = segextractAll(count:(count+originalinds(j)-1));
                    count = (count+originalinds(j));
                end
            end
            %%
            
        end
        
        
    end
    
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING

%% changing something in ALLBRANCH (e.g. location fo neurosn)
clear all; close all;
strname = 'xaa_Algn2Ons1_27Oct2017_1114_XLMAN25ms';

% =====================================================
load(['ALLBRANCHv2_' strname]);

bname = ALLBRANCH.SummaryStruct.birds(5).birdname;
assert(strcmp(bname, 'or74bk35'), 'asdfasd');

% --------- replace location
nneurons = length(ALLBRANCH.SummaryStruct.birds(5).neurons);
for i=1:nneurons
    
    ALLBRANCH.SummaryStruct.birds(5).neurons(i).NOTE_Location = 'LMAN';
end

% ------- SAVE
save(['ALLBRANCHv2_' strname], 'ALLBRANCH', '-v7.3');













