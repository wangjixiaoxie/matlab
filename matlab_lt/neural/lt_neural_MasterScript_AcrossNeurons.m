% %% neuron database [NONLEARNING (CODING)]
% clear NeuronDatabase;
% NeuronDatabase.global.basedir='/bluejay4/lucas/birds/bk7/NEURAL';
% ind=0;

%% TO DO LIST

% 1) FR - pitch correaltion, do for each time bin peri-syllable
% 2) 

%% ++++++++++++++++ OLD METHOD (NEURON DATABASE)
%% NOTE: CAN CURRENTLY ONLY DO ONE BIRD AT A TIME!!
clear all; close all;

%% ===== bk7

bk7_analysis_NeuronDatabase

%% ===== pk17

pk17gr57_analysis_NeuronDatabase


%% ====== bu77

bu77wh13_analysis_NeuronDatabase

%% ++++++++++++++++++++++++++++++++ NEW METHOD (SUMMARY STRUCT)

[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database();


BirdsToKeep = {'or74bk35'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {};
ExptToKeep = {};
RecordingDepth= [2490];
LearningOnly=0;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, BrainArea, ExptToKeep, ...
    RecordingDepth, LearningOnly);


%% ======== CONVERT TO TABLE
NeuronDatabase_table=struct2table(NeuronDatabase.neurons);
summary(NeuronDatabase_table);
disp(NeuronDatabase_table)

f = figure('Position',[200 200 1000 500]);
dat=table2cell(NeuronDatabase_table);
cnames=fieldnames(NeuronDatabase_table);
t = uitable('Parent', f, 'Data',dat,'ColumnName',cnames, 'Position',[0 0 900 450]);

% === if there are any empty neurons, then move another neuron


%% +++++++++++++++++++++ PLOTS
%% ===== PLOT A RANDOM SONG FILE (INCLUDING RAW NEURAL AND SPIKES)


%% ==== plot  spike waveforms for all neurons
close all;
numrandspikes=100;

lt_neural_MultNeur_SpkWaveforms(NeuronDatabase, numrandspikes);


%% FOR ALL NEURONS PLOT RASTERS, ETC, FOR A GIVEN MOTIF
if (0)
    % NOTE: USE BELOW INSTEAD!!!!!
    % need to throw this out (redundant with below) or put into a
    % function
    % ================ MOTIF STATISTICS (E.G. FIRING RATE, BURSTS, ...)
    for i=1:NumNeurons
        cd(NeuronDatabase.global.basedir);
        
        % - find day folder
        dirdate=NeuronDatabase.neurons(i).date;
        tmp=dir([dirdate '*']);
        assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
        cd(tmp(1).name);
        
        % - load data for this neuron
        batchf=NeuronDatabase.neurons(i).batchfile;
        channel_board=NeuronDatabase.neurons(i).chan;
        [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
        
        % --- EXTRACT DAT
        regexpr_str='WHOLEBOUTS';
        predur=6; % sec
        postdur=6; % sec
        alignByOnset=1;
        WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
        
        
        % ==================== Plot individually for this neuron
        useRescaled=0; % 1, then need to run LinTimeWarp first (plots scaled spikes, not song dat)
        plotAllSegs=0; % then plots each trial own plot.
        [Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs);
        
    end
end

%% ==== SONG MODULATION [SUMMARY OF FR AND ISI]
% ---- MORE ACCURATETLY, OUTPUTS FR AND ISI FOR USER DEFINED WINDOWS, FOR A
% GIVEN MOTIF
close all;
Window_relOnset={};
Window_relOffset={};

% === set params
regexpr_str='WHOLEBOUTS';
predur=6; % sec
postdur=6; % sec
alignByOnset=1;
WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps

% -- for summary analyses
motifwind=[-0.04 -0.04]; % e.g. 40ms premotif onset to 40ms pre motif offest
Window_relOnset{1}=[-6 -5.5]; % window rel 1st syl onset
Window_relOnset{2}=[-0.6 -0.1]; % window rel 1st syl onset
Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
Window_relOffset{2}=[5.5 6];

lt_neural_MultNeur_SongMod(NeuronDatabase, regexpr_str, predur, postdur, ...
    alignByOnset, WHOLEBOUTS_edgedur, Window_relOnset, motifwind, Window_relOffset)


%% COLLECT AND PLOT STATS ACROSS NEURONS [CHOOSE ONE MOTIF]
close all;
% ===================== OBSOLETE ---> USE MULTIPLE MOTIFS CODE BELOW
% 
% % ================ MOTIF STATISTICS (E.G. FIRING RATE, BURSTS, ...)
% motif_regexpr_str='[^b](b)bb';
% motif_regexpr_str='n(h)hh';
% motif_regexpr_str='nn(h)hhh[^h]';
% motif_regexpr_str='(h)hhhh-';
% 
% motif_regexpr_str='(g)b';
% motif_regexpr_str='(g)h';
% 
% motif_regexpr_str='(h)';
% 
% motif_regexpr_str='n(h)';
% motif_predur=0.2;
% motif_postdur=0.1;
% 
% % --- to plot entire bouts
% motif_regexpr_str='WHOLEBOUTS';
% motif_predur=6; % to be able to collect about 2 sec pre and post after stretch, account for up to 3x contraction (so get 6s)
% motif_postdur=6;
% % TO DO: make this not have to do linear stretch
% 
% lt_neural_MultNeur_MotifRasters(NeuronDatabase, motif_regexpr_str, motif_predur, motif_postdur)


%% SAME AS ABOVE, BUT CHOOSE MULTIPLE MOTIFS
% linearly warps across all motifs all neurons
close all;

% motif_regexpr_str={'(b)', '(h)'};
% motif_predur=0.2;
% motif_postdur=0;
% 
% motif_regexpr_str={'g(h)h'};
% motif_predur=0.1;
% motif_postdur=0.1;
% 
% 
% motif_regexpr_str={'n(h)hh', 'g(h)hh'};
% motif_predur=0.1;
% motif_postdur=0.1;
% 
% motif_regexpr_str={'[^vb](b)b', '[v](b)b'};
% motif_predur=0.2;
% motif_postdur=0.2;
% 
% motif_regexpr_str={'(y)'};
% motif_predur=0.1;
% motif_postdur=0.1;
 
motif_regexpr_str={'g(h)h', 'h(h)h', 'n(h)h'};
motif_predur=0.15;
motif_postdur=0.2;

% motif_regexpr_str={'r(g)', 'v(g)', 'h(g)', 'b(g)', 'n(g)'};
% motif_predur=0.05;
% motif_postdur=0.05;

% motif_regexpr_str={'r(g)', 'v(g)', 'h(g)', 'b(g)', 'n(g)'};
% motif_predur=0.05;
% motif_postdur=0.05;

% motif_regexpr_str={'(b)[sr]', '(b)g'}; % USEFUL COMPARISON?
% motif_predur=0.2;
% motif_postdur=0.1;

% motif_regexpr_str={'(k)e', '(k)k', '(k)l'}; % 
% motif_predur=0.1;
% motif_postdur=0.05;

% motif_regexpr_str={'nn(h)hh', 'nnh(h)h'};
% motif_predur=0.3;
% motif_postdur=0.2;
% 
% motif_regexpr_str={'g(h)', 'n(h)', 'g(v)', 'g(b)'};
% motif_predur=0.1;
% motif_postdur=0.1;
% 
% ------------------------------------ pk17
% motif_regexpr_str={'g(h)', 'n(h)', 'g(v)', 'g(b)'};
% motif_regexpr_str={'k(d)ccbgh', 'w(d)ccbgh'};
% motif_predur=0.2;
% motif_postdur=0.1;

% motif_regexpr_str={'g(b)', 'c(b)'};
% motif_predur=0.2;
% motif_postdur=0.1;


% ------------------------- BU77
% motif_regexpr_str={'a(b)bh', 'j(b)h'};
% motif_predur=0.2;
% motif_postdur=0.3;

% ------------------ WH6
% motif_regexpr_str={'nl(c)chbg', 'ajkl(c)chbg', 'amksd(v)b'};
% motif_regexpr_str={'nl(c)chbg', 'ajkl(c)chbg', 'amksdv(b)', 'cch(b)', 'd(m)', 'a(m)'};
% motif_regexpr_str={'nl(c)c', 'jkl(c)c', 'nlc(c)', 'jklc(c)'};
% motif_regexpr_str={'v(b)', 'klcch(b)', 'nlcch(b)'};
% % motif_regexpr_str={'d(m)', 'a(m)'};
% motif_predur=0.2;
% motif_postdur=0.1;

% ------ BR92
motif_regexpr_str={'nk(h)', 'dd(d)', 'ddd(h)', 'ag(c)c'};
motif_regexpr_str={'nk(h)', 'd(d)d'};
motif_regexpr_str={'nk(h)', 'ag(c)'};
% motif_regexpr_str={'h(d)', '-(d)'};

motif_predur=0.2;
motif_postdur=0.1;

% ----- OR74
motif_regexpr_str={'ab(g)', 'anb(g)'};
motif_predur=0.35;
motif_postdur=0.1;



LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)
WHOLEBOUTS_edgedur = '';

% --- for whole motif
% motif_regexpr_str={'WHOLEBOUTS'};
% motif_predur=6; % to be able to collect about 2 sec pre and post after stretch, account for up to 3x contraction (so get 6s)
% motif_postdur=6;
% LinScaleGlobal=1;
% WHOLEBOUTS_edgedur = 6;
% --

lt_neural_MultNeur_MotifRasters_v2(NeuronDatabase, motif_regexpr_str, motif_predur, motif_postdur, LinScaleGlobal, WHOLEBOUTS_edgedur);


%% ======= PLOT ALL SYLS IN MOTIFS (SAME AS ABOVE, BUT AUTOMATICALLY PICKS OUT SYLS)


%% ======== ARRANGE RESPONSES BASED ON CHANNEL AND DEPTH
% PLOT RESPONSE arranged by channel, depth, and date.
close all;

motif_regexpr_str={'j(b)'};
motif_predur=0.2;
motif_postdur=0.1;

LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)
WHOLEBOUTS_edgedur = '';

lt_neural_MultNeur_ArrangeResponses(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, WHOLEBOUTS_edgedur);


%% ======== EXTRACT ALL ONE-BACK TRANSITIONS (including probabilities and frequencies, both div and conv)
close all;

[TRANSMATRIX] = lt_neural_MultNeur_GetTrans(NeuronDatabase);


%% ======== MATRIX OF MEAN RESPONSE FOR ALL ONE-BACK TRANSITIONS
% --- i.e. for each neuron, plots all conv and div branch points and
% smoothed firing rates

close all;
mindatFreq=5; % only analyze transitions with N or more cases
lt_neural_MultNeur_AllOnebackNeur(NeuronDatabase, TRANSMATRIX, mindatFreq)


%% ========= FOR EACH BRANCH POINT, PLOT RESPONSE BY EACH NEURON - 
% ALSO PLOT "PSEUDOPOPULATION RESPONSE"
% ALSO PLOT TIMECOURSE OF DEVIATION ACROSS TRANSITIONS
% KEY THINGS ABOUT THIS ANAYLSYSI:"
% % ============= SINGLE BIRD ANALYSIS, assumes all neurons have same
% sylalbles.


close all;
branchtype = 'conv';
plotOnlySummary = 1 ; %then will auto skip or close figs other than summaries.
lt_neural_MultNeur_pseudopop(NeuronDatabase, TRANSMATRIX, branchtype, plotOnlySummary);


%% ====== REVISED VERSION OF ABOVE - CAN COMPARE ANY BRANCHES, NOT JUST ONE BACK
% KEY ADVANTAGE: DOES NEG AND POSITIVE CONTROLS, MATCHING SAMPLE SIZES

close all;
closeEarlyFigs = 1;
saveFigs = 0; % NOTE: very large files!!! keep off.
NCycles = 10; % cycles to perform shuffles
lt_neural_MultNeur_AnovaTcourse(NeuronDatabase, NCycles, closeEarlyFigs)

% TO DO: 
% each time point, compare cdf of anova.
% each time point, plot pdf of pvalues (shuffle comparison)


%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ========================== BUILD LINEAR MODEL OF FIRING RATES
close all;

% ========================= FREQ WINDOWS 
% 1) bk7
FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
            'v', [2450 4300]};
FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
            'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)

% use below for WN on h
FFparams.cell_of_FFtimebins={'h', [0.034 0.038], 'b', [0.053 0.07], ...
            'v', [0.052 0.07]}; % WN on g H


% 2) bu77
FFparams.cell_of_freqwinds={'b', [2700 3900], 'h', [2600 3900], 'a', [1300 2600]};
FFparams.cell_of_FFtimebins={'b', [0.031 0.04], 'h', [0.038 0.052], 'a', [0.056 0.081]}; % in sec, relative to onset (i.e. see vector T)


% 3) wh6
FFparams.cell_of_freqwinds={'c', [2100 3100], 'h', [2800 4000], 'b', [2700 3800], ...
    'a', [1300 2200], 's', [4000 5100], 'd', [900 2000],  'n', [3300 4300], 'v', [2600 4000]};
FFparams.cell_of_FFtimebins={'c', [0.035 0.05], 'h', [0.023 0.033], 'b', [0.034 0.038], ...
    'a', [0.06 0.08], 's', [0.02 0.023], 'd', [0.02 0.027],  'n', [0.029 0.05], 'v', [0.028 0.036]}; % in sec, relative to onset (i.e. see vector T)

% 4) br92 - TEMPORARY
FFparams.cell_of_freqwinds={'a', [800 1400], 'c', [1300 1800], 'h', [2300 3600], 'd', [1300 3400], 'k', [800 1800]};
FFparams.cell_of_FFtimebins={'a', [0.04 0.049], 'c', [0.04 0.049], 'h', [0.040 0.049], 'd', [0.032 0.052], 'k', [0.05 0.055]}; % in sec, relative to onset (i.e. see vector T)

saveOn = 0; % then saves SYLNEURDAT
plotSpec = 1; % to plot raw spec overlayed with PC and windows.
plotOnSong = 15; % will only start plotting spec once hit this song num.
plotSyl = ''; % to focus on just one syl. NOT DONE YET

SYLNEURDAT = lt_neural_MultNeur_CollectFeats(NeuronDatabase, FFparams, saveOn, ...
    plotSpec, plotOnSong, plotSyl);

% NOTE DOWN PARAMS OF LAST RUN:
% good (all neurons)

%% ================= MOVING WINDOW CORRELATION BETWEEN FR AND FF




%% ======= LOAD LAST SAVED STRUCT

savedir='/bluejay4/lucas/birds/bk7/SUMMARYDAT_MultNeur_CollectFeats';
cd(savedir)
load('SYLNEURDAT')


%% ========== construct model 
close all;

lt_neural_MultNeur_MixedEffects(SYLNEURDAT);


%% +++++++++++++++++++++++++++++++++++++++++++++++++ NEWER STUFF



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++ LEARNING
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% IMPORTANT - THIS CODE IS FOR ONE NEURON ONLY - WILL TAKE FIRST NEURON IN NEURONDATABASE

% ==========================================================BK7
% FOR EACH NEURON, PLOT RASTER AND SMOOTHED FIRING FOR A GIVEN MOTIF ACROSS
% TIME allsyls
% - I.E. same as above code, except noting when learning began, gets FF,
% hit/escape, and catch song information
close all; 
% motif_regexpr_str={'g(h)'};
motif_regexpr_str={'g(h)', 'nnn(h)h'};
motif_regexpr_str={'g(h)', 'nn(h)h','v(b)bb', 'g(b)'};
% motif_regexpr_str={'g(h)', 'nn(h)h','v(b)b', 'n(b)b', 'jk(k)kkkl', 'g(v)', 'g([bv])', '[oh](p)pp', 'r(s)'};
% motif_regexpr_str={'g(h)', 'nn(h)h','v(b)b', 'g([vb])', 'j(k)kk', 'y(o)o', 'o(p)p'}; % latest
% motif_regexpr_str={'g(h)', 'nnn(h)hh','v(b)bb', 'n(b)bb', 'jk(k)kkkl', 'g(v)', 'o(p)pp', 'h(p)pp', 'r(s)'}; % -- old

motif_predur=0.1;
motif_postdur=0.1;
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

% NOTE: for below, run linear model stuff above
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will 
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
    % +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
    % not required (can set as []);
% NOTE: will also determine whether was hit or miss, based on WN sound
% detection.
        
% LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'}; % LearnLMAN1 - first epoch
% WNchangeDateStrings={'06Oct2016-1225'}; % LearnLMAN1 - second epoch
% WNchangeDateStrings={'17Oct2016-1332'}; % LearnLMAN2
% WNchangeDateStrings={'17Oct2016-1940'}; % LearnLMAN2 - second epoch - arbitrary, since did not actually change WN during this.


OnlyPlotNoHit=0; % then only plots trials that were not hit (WN)
TrialBinSize=10;


% ================================================================== bu77

% FOR EACH NEURON, PLOT RASTER AND SMOOTHED FIRING FOR A GIVEN MOTIF ACROSS
% TIME allsyls
% - I.E. same as above code, except noting when learning began, gets FF,
% hit/escape, and catch song information
close all; 
motif_regexpr_str={'a(b)', 'j(b)'};
motif_regexpr_str={'a(b)', 'j(b)', 'jb(h)', 'ab(b)', 'abb(h)', '(a)b', '(j)b'};
motif_regexpr_str={'a(b)'};

motif_predur=0.4;
motif_postdur=0.2; 
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

% NOTE: for below, run linear model stuff above
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will 
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
    % +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
    % not required (can set as []);
% NOTE: will also determine whether was hit or miss, based on WN sound
% detection.
        
% LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'}; % LearnLMAN1 - first epoch
% WNchangeDateStrings={'06Oct2016-1225'}; % LearnLMAN1 - second epoch
% WNchangeDateStrings={'17Oct2016-1332'}; % LearnLMAN2
% WNchangeDateStrings={'17Oct2016-1940'}; % LearnLMAN2 - second epoch - arbitrary, since did not actually change WN during this.


OnlyPlotNoHit=0; % then only plots trials that were not hit (WN)
TrialBinSize=15; % default 10

% ================================================================== wh6

% FOR EACH NEURON, PLOT RASTER AND SMOOTHED FIRING FOR A GIVEN MOTIF ACROSS
% TIME allsyls
% - I.E. same as above code, except noting when learning began, gets FF,
% hit/escape, and catch song information
close all; 
motif_regexpr_str={'nl(c)chbg', 'ajkl(c)chbg', 'amksd(v)b'}; % LMANlearn1
motif_regexpr_str={'nlcch(b)', 'klcch(b)', 'dv(b)', 'l(c)', 'lc(c)', 'lcc(h)', ...
    'k(s)', 'd(v)', 'g(a)'}; % wh6
motif_regexpr_str={'nlcch(b)g', 'klcch(b)g', 'amksdv(b)'}; % LMANlearn2

motif_predur=0.5;
motif_postdur=0.2; 
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

% NOTE: for below, run linear model stuff above
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will 
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
    % +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
    % not required (can set as []);
% NOTE: will also determine whether was hit or miss, based on WN sound
% detection.
        
% LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'}; % LearnLMAN1 - first epoch
% WNchangeDateStrings={'06Oct2016-1225'}; % LearnLMAN1 - second epoch
% WNchangeDateStrings={'17Oct2016-1332'}; % LearnLMAN2
% WNchangeDateStrings={'17Oct2016-1940'}; % LearnLMAN2 - second epoch - arbitrary, since did not actually change WN during this.


OnlyPlotNoHit=0; % then only plots trials that were not hit (WN)
TrialBinSize=20; % default 10



% === ONLY DOES ONE NEURON AT A TIME!! - FIRST INDEX IN NeuronDatabase;
lt_neural_MultNeur_MotifRasters_Learning(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, OnlyPlotNoHit, TrialBinSize)

% IN PROGRESS!!:: LinScaleGlobal.

% ================================================================== wh6

% FOR EACH NEURON, PLOT RASTER AND SMOOTHED FIRING FOR A GIVEN MOTIF ACROSS
% TIME allsyls
% - I.E. same as above code, except noting when learning began, gets FF,
% hit/escape, and catch song information
close all; 
motif_regexpr_str={'ddd(h)', 'nk(h)'};

motif_predur=0.4;
motif_postdur=0.1; 
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

% NOTE: for below, run linear model stuff above
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will 
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
    % +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
    % not required (can set as []);
% NOTE: will also determine whether was hit or miss, based on WN sound
% detection.
        
% LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'}; % LearnLMAN1 - first epoch
% WNchangeDateStrings={'06Oct2016-1225'}; % LearnLMAN1 - second epoch
% WNchangeDateStrings={'17Oct2016-1332'}; % LearnLMAN2
% WNchangeDateStrings={'17Oct2016-1940'}; % LearnLMAN2 - second epoch - arbitrary, since did not actually change WN during this.


OnlyPlotNoHit=0; % then only plots trials that were not hit (WN)
TrialBinSize=20; % default 10



% === ONLY DOES ONE NEURON AT A TIME!! - FIRST INDEX IN NeuronDatabase;
lt_neural_MultNeur_MotifRasters_Learning(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, OnlyPlotNoHit, TrialBinSize)

% IN PROGRESS!!:: LinScaleGlobal.
%% ============== LEARNING CONTOURS SUMMARY

% === GIVEN ACTUAL MOTIFS FOR THIS BIRD, DETERMINE WHICH SYLS TO LOOK AT
MotifsActual = {'nlcchb', 'jklcchb', 'ga', 'mksd'};
motif_regexpr_str = {}; % will put motifs for extraction here

for i=1:length(MotifsActual)
   motif_actual = MotifsActual{i};
   
   % for each vocalization in the motif, extract one segment
   numvocals = length(motif_actual);
   
   for ii=1:numvocals
       if ii==1
       segmentmotif = ['(' motif_actual(ii) ')' motif_actual(ii+1:end)];
          
       elseif ii == numvocals
       segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')'];
                     
       else
           
       segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')' motif_actual(ii+1:end)];
           
       end
       
%        disp(segmentmotif);

motif_regexpr_str = [motif_regexpr_str segmentmotif];        
        
   end
end

% ===============
close all; 
% bk7
motif_regexpr_str={'g(h)', 'nn(h)h','v(b)bb', 'g(b)'};

% - wh6
motif_regexpr_str={'nl(c)chbg', 'ajkl(c)chbg', 'amksd(v)b'}; % LMANlearn1
motif_regexpr_str={'nlcch(b)', 'klcch(b)', 'dv(b)', 'l(c)', 'lc(c)', 'lcc(h)', ...
    'k(s)', 'd(v)', 'g(a)'}; % wh6
motif_regexpr_str={'nlcch(b)g', 'klcch(b)g', 'amksdv(b)'}; % LMANlearn2
motif_regexpr_str={'nlcch(b)', 'klcch(b)', 'dv(b)', 'l(c)', 'lc(c)', 'lcc(h)', ...
    'k(s)', 'd(v)', 'g(a)', 'm(n)'}; % wh6

% br92br54
motif_regexpr_str={'ddd(h)', 'nk(h)'};
motif_regexpr_str={'c(d)dd', 'h(d)dd', 'ag(c)c', 'n(k)h', 's(d)'};


% or74
motif_regexpr_str={'a(b)g', 'an(b)'};
motif_regexpr_str={'ab(g)', 'anb(g)'};

% bu77
motif_regexpr_str={'a(b)g', 'an(b)'};
motif_regexpr_str={'a(b)', '(a)b', 'j(b)'};
motif_regexpr_str={'a(b)'};


motif_predur=0.2;
motif_postdur=0.1; 
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

% NOTE: for below, run linear model stuff above
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will 
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
    % +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
    % not required (can set as []);
% NOTE: will also determine whether was hit or miss, based on WN sound
% detection.
        
% LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'}; % LearnLMAN1 - first epoch
% WNchangeDateStrings={'06Oct2016-1225'}; % LearnLMAN1 - second epoch
% WNchangeDateStrings={'17Oct2016-1332'}; % LearnLMAN2
% WNchangeDateStrings={'17Oct2016-1940'}; % LearnLMAN2 - second epoch - arbitrary, since did not actually change WN during this.


OnlyPlotNoHit=0; % then only plots trials that were not hit (WN)
TrialBinSize=20; % default 10
DivisorBaseSongs = 1;


% === ONLY DOES ONE NEURON AT A TIME!! - FIRST INDEX IN NeuronDatabase;
lt_neural_MultNeur_MotifRasters_Learning2(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, OnlyPlotNoHit, ...\
    TrialBinSize, DivisorBaseSongs)

% IN PROGRESS!!:: LinScaleGlobal.



%% ====== SUMMARY OF D-PRIME STUFF (BIN OVER TIME AND TRIALS) only plots summary figures
% d-prime of firing rate change in premotor window summary plots, relative
% to learning timepoints and FF change
close all;
% === NEW PARAMS.
% motif_regexpr_str={'a(b)', 'j(b)', 'jb(h)', 'ab(b)', 'abb(h)', '(a)b',
% '(j)b'}; % which bird
motif_regexpr_str={'nlcch(b)', 'klcch(b)', 'dv(b)', 'l(c)', 'lc(c)', 'lcc(h)', ...
    'k(s)', 'd(v)', 'g(a)'}; % wh6
UseEntireBaseline = 0; % if 1, uses entire baseline, otherwise uses 1st bin.
TypesToPlot={'mean'}; % for d-prime summary metric (which ones to plot?)
% TypesToPlot={'mean', 'abs', 'std', 'corr'}; % for d-prime summary metric (which ones to plot?)
LinScaleGlobal = 0; % NOTE: IN PROGRESS
TrialBinSize=10;
premotor_wind=[-0.1 0.04]; % in sec, relative to onset of token in motif.

lt_neural_MultNeur_MotifRasters_LearnSum(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, ...
    OnlyPlotNoHit, UseEntireBaseline, TypesToPlot, TrialBinSize, premotor_wind)


%% ============ SUMMARY, TRIAL BY TRIAL [IMPORTANT]

close all;

% motif_regexpr_str={'a(b)', 'j(b)', 'jb(h)', 'ab(b)', 'abb(h)', '(a)b', '(j)b'};
% motif_regexpr_str={'a(b)'};
% motif_regexpr_str={'a(b)', 'j(b)'};
motif_regexpr_str={'jkl(c)', 'nl(c)', 'cc(h)', 'cch(b)', 'ksd(v)', 'sdv(b)', 'm(n)', 'ks(d)', 'g(a)'}; % bu77?
motif_regexpr_str={'nlcch(b)', 'klcch(b)', 'dv(b)', 'l(c)', 'lc(c)', 'lcc(h)', ...
    'k(s)', 'd(v)', 'g(a)'}; % wh6
motif_predur=0.4;
motif_postdur=0.2; 

LinScaleGlobal = 0; % NOTE: IN PROGRESS
premotor_wind=[-0.1 0.03]; % in sec, relative to onset of token in motif.

lt_neural_MultNeur_LearningAnaly(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, premotor_wind);



%% ====== ACUTE EFFECT OF WN?
% NOTE: RUN ABOVE PARAMS FIRST
close all;
motif_regexpr_str={'g(h)h'}; % bk7
motif_regexpr_str={'a(b)bh'}; % bu77
motif_regexpr_str={'lcch(b)g'}; % bu77
motif_regexpr_str={'nk(h)'}; % br92
motif_regexpr_str={'h(d)dd', 'c(d)dd'}; % br92
motif_regexpr_str={'ab(g)'}; % or74

motif_predur=0.2;
motif_postdur=0.1;
NumTrialsBin=20;
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

lt_neural_MultNeur_MotifRasters_WNacute(NeuronDatabase, motif_regexpr_str, ...
    motif_predur,motif_postdur, FFparams, NumTrialsBin, LinScaleGlobal)







%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ RANDOM THINGS

%% remove song dat from metadat
% THIS IS UP TO DATE. WILL SAVE SONG IN UPPER DIR AND REMOVE SONG DAT FROM
% META DAT. RUNNING WILL NOT CORRUPT ANYTHING. 
% RUN THIS FOR NEW FINALIZED NEURONS.

numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        disp(['========= bird ' num2str(i) ' neuron ' num2str(ii) ', ' ...
            SummaryStruct.birds(i).neurons(ii).batchfilename ' (if blank did nothing)']);
        
        Datstruct = SummaryStruct.birds(i).neurons(ii);
        cd(Datstruct.dirname)
        
        % 1) === IS THIS BATCH'S SONG ALREADY SAVED UP ONE LEVEL?
        cd ..
        songfname = ['SongDat_' Datstruct.batchfilename '.mat'];
        songAlreadyExists = 0;
        if exist(songfname, 'file')==2
            % then already exists, skip
            songAlreadyExists = 1;
        end
        cd(Datstruct.dirname)
        
        % 2) ====  Extract and save song if necessary.
        if songAlreadyExists == 0
            disp('EXTRACTING SONG ...');
            % -- if not saved, then save
            metadat = load('MetaDat');
            
            % save song data in cell array
            numsongs = length(metadat.metaDat);
            SongCellArray = {};
            for j=1:numsongs
                SongCellArray = [SongCellArray single(metadat.metaDat(j).songDat)];
            end
            
            % save
            cd ..
            save(songfname, 'SongCellArray', '-v7.3');
            disp(['-- extracted and saved ! (' songfname ')']);
            cd(Datstruct.dirname)            
        end
        
        % 3) ===== if Song not removed from MetaDat, then do that.
        if exist('DONE_RemovedSongDatFromMetaDat', 'file') ==2
            % then already removed, DO NOTHING.
            
        else
            disp('REMOVING SONG FROM METADAT');
            metadat = load('MetaDat');
            % ------ 2) don't yet remove from metaDat (for backwards compatibility testing)
            % INSTEAD CHANGE THE NAME. IF WORKS FINE, THEN REMOVE.
            metaDat = rmfield(metadat.metaDat, 'songDat'); % also name change
            
            % move old metadat to new file name
            % DELETE OLD METADAT
            eval('!rm MetaDat.mat');
            disp('DELETED OLD METADAT!');
            
            % save new
            save('MetaDat.mat', 'metaDat');
            
            fid = fopen('DONE_RemovedSongDatFromMetaDat', 'w');
            fclose(fid);
        end
    end
end







