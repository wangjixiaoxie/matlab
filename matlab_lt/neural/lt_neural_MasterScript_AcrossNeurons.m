% %% neuron database [NONLEARNING (CODING)]
% clear NeuronDatabase;
% NeuronDatabase.global.basedir='/bluejay4/lucas/birds/bk7/NEURAL';
% ind=0;

%% TO DO LIST

% 1) FR - pitch correaltion, do for each time bin peri-syllable
% 2) 

%% NOTE: CAN CURRENTLY ONLY DO ONE BIRD AT A TIME!!
clear all; close all;

%% ===== bk7

bk7_analysis_NeuronDatabase

%% ===== pk17

pk17gr57_analysis_NeuronDatabase


%% ====== bu77

bu77wh13_analysis_NeuronDatabase

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
NumNeurons=length(NeuronDatabase.neurons);

figcount=1;
subplotrows=NumNeurons;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

plotcols=lt_make_plot_colors(NumNeurons, 0, 0);

for i=1:NumNeurons
    
    cd(NeuronDatabase.global.basedir);
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    
    inds=find(NeurDat.spikes_cat.cluster_class(:,1)==NeuronDatabase.neurons(i).clustnum); % get desired clust
    % get random subsamp
    inds=inds(randperm(length(inds), numrandspikes)); % get subset of spikes
    spkwaves=NeurDat.spikes_cat.spikes(inds, :);
    fs=NeurDat.metaDat(1).fs;
    
    % -- plot individual spikes
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['neuron ' num2str(i)]);
    tt=1:size(spkwaves,2);
    tt=1000*tt/fs;
    plot(tt, spkwaves', 'Color', plotcols{i});
    
    % -- overlay mean + std
    spkwaves_mean=mean(spkwaves,1);
    spkwaves_std=std(spkwaves,0, 1);
    plot(tt, spkwaves_mean, 'Color', plotcols{i}, 'LineWidth', 3);
    
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        % then is mu
        lt_plot_text(tt(1), 200, 'MU');
        
    end
    
    axis tight
end


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


%% ==== MOTIF MODULATION

% ========== FOR A GIVEN NEURON, COMPARE RASTERS FOR TWO MOTIFS, ALIGNING
% THOSE MOTIFS


% ============= GIVEN A MOTIF, ALIGN RASTERS FOR ALL NEURONS (that have
% data for that motif)

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
 
% motif_regexpr_str={'g(h)', 'h(h)', 'n(h)'};
% motif_predur=0.1;
% motif_postdur=0.1;

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
motif_regexpr_str={'abb(h)', 'jb(h)'};
motif_predur=0.2;
motif_postdur=0.3;
% 
% 

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

%% ======== ARRANGE RESPONSES BASED ON CHANNEL AND DEPTH
% PLOT RESPONSE arranged by channel, depth, and date.

motif_regexpr_str={'a(b)', 'j(b)'};
motif_predur=0.1;
motif_postdur=0.2;

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
lt_neural_MultNeur_AllOnebackNeur(NeuronDatabase, TRANSMATRIX)


%% ========= FOR EACH BRANCH POINT, PLOT RESPONSE BY EACH NEURON - 
% ALSO PLOT "PSEUDOPOPULATION RESPONSE"
% ALSO PLOT TIMECOURSE OF DEVIATION ACROSS TRANSITIONS

close all;
branchtype = 'conv';
plotOnlySummary = 1 ; %then will auto skip or close figs other than summaries.
lt_neural_MultNeur_pseudopop(NeuronDatabase, TRANSMATRIX, branchtype, plotOnlySummary);


%% ====== REVISED VERSION OF ABOVE - CAN COMPARE ANY BRANCHES, NOT JUST ONE BACK
close all;
closeEarlyFigs = 1;
NCycles = 1000; % cycles to perform shuffles
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
% FFparams.cell_of_FFtimebins={'h', [0.034 0.038], 'b', [0.053 0.07], ...
%             'v', [0.052 0.07]}; % WN on g H


% 2) bu77
FFparams.cell_of_freqwinds={'b', [2700 3900], 'h', [2600 3900], 'a', [1300 2600]};
FFparams.cell_of_FFtimebins={'b', [0.031 0.04], 'h', [0.038 0.052], 'a', [0.056 0.081]}; % in sec, relative to onset (i.e. see vector T)


saveOn = 1; % then saves SYLNEURDAT
plotSpec = 1; % to plot raw spec overlayed with PC and windows.
SYLNEURDAT = lt_neural_MultNeur_CollectFeats(NeuronDatabase, FFparams, saveOn, plotSpec);

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


% === ONLY DOES ONE NEURON AT A TIME!! - FIRST INDEX IN NeuronDatabase;
lt_neural_MultNeur_MotifRasters_Learning(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, OnlyPlotNoHit, TrialBinSize)

% IN PROGRESS!!:: LinScaleGlobal.

%% ====== SUMMARY OF D-PRIME STUFF (BIN OVER TIME AND TRIALS) only plots summary figures
% d-prime of firing rate change in premotor window summary plots, relative
% to learning timepoints and FF change
close all;
% === NEW PARAMS.
motif_regexpr_str={'a(b)', 'j(b)', 'jb(h)', 'ab(b)', 'abb(h)', '(a)b', '(j)b'};
% motif_regexpr_str={'a(b)'};
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

motif_regexpr_str={'a(b)', 'j(b)', 'jb(h)', 'ab(b)', 'abb(h)', '(a)b', '(j)b'};
motif_regexpr_str={'a(b)'};
% motif_regexpr_str={'a(b)', 'j(b)'};
motif_predur=0.15;
motif_postdur=0.2; 

LinScaleGlobal = 0; % NOTE: IN PROGRESS
premotor_wind=[-0.03 0.01]; % in sec, relative to onset of token in motif.

lt_neural_MultNeur_LearningAnaly(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, FFparams, premotor_wind);



%% ====== ACUTE EFFECT OF WN?
% NOTE: RUN ABOVE PARAMS FIRST
close all;
motif_regexpr_str={'g(h)h'}; % bk7
motif_regexpr_str={'a(b)bh'}; % bu77
motif_predur=0.25;
motif_postdur=0.4;
NumTrialsBin=20;
LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

lt_neural_MultNeur_MotifRasters_WNacute(NeuronDatabase, motif_regexpr_str, ...
    motif_predur,motif_postdur, FFparams, NumTrialsBin, LinScaleGlobal)







