%% lt 11/17/17 - to analyze smoothed firing across multiple channels
%% ##################################### 1) EXTRACT DATA
%% ====================== SWITCH ONE (11:25a)
clear all; close all;
% ----- different datasets, each with different batch/dir combination
basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110517_RALMANlearn2/';
subdirs = {'', ''}; % these must be dirs within basedir. they will be used as fieldnames in extracted structure
Batchfiles = {'BatchSw1Pre', 'BatchSw1Post'}; % must be one for each subdir;
Conditions = {'base', 'WNon'};

% ------ params for all subdirs
chanstoplot = [9 14 17 21]; % chip channels. leave empty to get all
motifstoplot = {'(a)ab', 'a(a)b', 'a(b)', '(j)jb', 'j(j)b', 'jj(b)', ...
    'jb(h)', 'jbh(h)', 'h(g)'}; % cell arary of motifs, strings
% sylstoalign = [2 2 3 4 1]; % which syl in that motif? (one for each motif
pretime = 0.1; % sec, from onset
posttime = 0.1; % sec, from offset
plotRaw = 0; % plot each extracted trial (raw, smoothed, and spec)


DATSTRUCT = lt_neural_BatchSmthEXTRACT(basedir, subdirs, Batchfiles, chanstoplot, motifstoplot, ...
    '', pretime, posttime, plotRaw, Conditions);



%% ====================== SWITCH TWO (a)
close all;
% ----- different datasets, each with different batch/dir combination
basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110517_RALMANlearn2/';
subdirs = {'', ''}; % these must be dirs within basedir. they will be used as fieldnames in extracted structure
Batchfiles = {'BatchSw2Pre', 'BatchSw2Post'}; % must be one for each subdir;
Conditions = {'WNon', 'WNoff'};


DATSTRUCT = lt_neural_BatchSmthEXTRACT(basedir, subdirs, Batchfiles, chanstoplot, motifstoplot, ...
    '', pretime, posttime, plotRaw, Conditions);



%% ######################################## CLEAN ALL SONG FILES FOR BATCHES
close all
% basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111817_Afternoon_DirUndir/';
% subdirs = {'DIR', 'UNDIR'}; % these must be dirs within basedir. they will be used as fieldnames in extracted structure
% Batchfiles = {'batchall', 'batchall'}; % must be one for each subdir;

basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110517_RALMANlearn2/';
subdirs = {'', ''}; % these must be dirs within basedir. they will be used as fieldnames in extracted structure
Batchfiles = {'BatchSw2Pre', 'BatchSw2Post'}; % must be one for each subdir;

pretime = 0.1; % sec, from onset
posttime = 0.1; % sec, from offset
chanstoplot = [14 21]; % chip channels. leave empty to get all
skipifdone =1; % if 1, then skips if already done and old notmat labels are identical to current and chans to plot identical. if 0, then always redos
plotspecgram =1;
freqrange = [1800 4000]; % plots mean power over this range. if 0, then doesn't plot

lt_neural_CleanSylDat(basedir, subdirs, Batchfiles, pretime, posttime, ...
    chanstoplot, skipifdone, plotspecgram, freqrange);

%% ###################################### PLOTS [FOR INDIVIDUAL SWITCHES]

% ======================== PLOT RAW
% go to:
lt_neural_BatchSmth_ALLPLOTS;

% ======================== XCOV OF FIRING BETWEEN CHANNELS
close all;
windowmax = 0.05;
binsize = 0.005;
premotorwind = [-0.035 0.025]; % use for ch14-21 on 11/12, morning.
% premotorwind = [-0.8 0.02];
chan1 = 14;
chan2 = 17;
mm = 1;
dirfield = 'WNon';

lt_neural_BatchSmthXCOV(DATSTRUCT, windowmax, binsize, premotorwind, chan1, ...
    chan2, mm, dirfield, pretime);


% ========================= PLOT
close all;
fs = 30000;
fnames  = fieldnames(DATSTRUCT);
chanstoplot = DATSTRUCT.(fnames{1}).params.chanstoplot;
motifstoplot = DATSTRUCT.(fnames{1}).params.Motifstoplot;
motifstoplot = {'a(b)', 'j(b)'};
lt_neural_BatchSmthPLOT(DATSTRUCT, chanstoplot, motifstoplot, ...
    '', pretime, posttime, plotRaw, fs)

%% ##################################### 2) COMBINE MULTIPLE SWITCHES BY HAND
clear all; close all;

basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/';
ListOfDatstructs = {...
    'DATSTRUCT_BatchSm_17Nov2017_1235.mat',...
    'DATSTRUCT_BatchSm_17Nov2017_1227.mat'};


DATAllSwitches = lt_neural_BatchSmthCOMPILE(basedir, ListOfDatstructs);


%% ################################# PLOTS [COMBINED MULT SWITCHES]
%% ====================== PLOT RAW FOR EACH CHANNEL
close all;
motiftoplot = 'j(b)';
chanstoplot = find(~cellfun('isempty', DATAllSwitches.switch(1).motif(1).batchinorder(1).DatAll));

for cc = chanstoplot
    for i=1:length(DATAllSwitches.switch)
        
        nummotifs = length(DATAllSwitches.switch(i).motif);
        
        for ii=1:nummotifs
            
            numbatches = length(DATAllSwitches.switch(i).motif(ii).batchinorder);
            figcount=1;
            subplotrows=6;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            if ~strcmp(DATAllSwitches.switch(i).motif(ii).batchinorder(1).motifname, ...
                    motiftoplot)
                continue
            end
            
            for bb = 1:numbatches
                
                
                motif = DATAllSwitches.switch(i).motif(ii).batchinorder(bb).motifname;
                cond = DATAllSwitches.switch(i).motif(ii).batchinorder(bb).condition;
                % ========= plot
                datmat = DATAllSwitches.switch(i).motif(ii).batchinorder(bb).DatAllRaw{cc};
                t = DATAllSwitches.switch(i).motif(ii).batchinorder(bb).t;
                
                for j = 1:size(datmat,1);
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    hsplots = [hsplots hsplot];
                    title(['ch' num2str(cc) ',' motif ',sw' num2str(i) ',[' cond ']']);
                    
                    plot(t, datmat(j,:), 'k');
                    axis tight;
                    ylim([-175 175]);
                end
                
            end
            
            %             linkaxes(hsplots, 'xy');
            pause;
            close all;
            
        end
        
    end
end

%% ============================== NOTE DOWN WHICH ARE NOISY TRIALS
% [NOTE: INSTEAD OF THIS CHECK NOISE ON ORIGINAL SONG FILES]
% ---- will plot raw data for each channel/motif, all trials. type the
% trials that are noise. will save a file noting down noise trials dor each
% channel/motif.
% ---- NOTE: the saved file will be over the original song files, so that this
% code does not have to be rerun for every dat struct

close all;
lt_neural_BatchSmth_Clean(DATAllSwitches)

%% ==== 1) PLOT, each channel, plot each trial + mean, comparing conditions
close all;
plotRaw = 1; % if 1, then plots all trials and blocks overlaied. if 0 then goes straight to cross corr summary.
motifstoplot = {}; % if empty, plots all
premotor_wind = [-0.03 0.02]; % for cross correlation (rel syl onset);
removeNoiseTrials = 0;

lt_neural_BatchSmth_Premotor(DATAllSwitches, motifstoplot, premotor_wind, ...
    plotRaw, removeNoiseTrials)




%% ==== 1) [CONTEXTUAL SEPARATION] MEAN FR, COMPARING SYLS, PRE AND POST
close all;

clear MotifSets;
% MotifSets{1} = {'a(b)', 'j(b)'};
% MotifSets{2} = {'ab(h)', 'jb(h)'};
MotifSets{1} = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)'};
% MotifSets{2} = {'(a)ab', 'a(a)b'};
% MotifSets{3} = {'(j)jb', 'j(j)b'};
% MotifSets{4} = {'jb(h)', 'jbh(h)'};

plotRaw =1; % if 1, then plots FR traces (USEFUL). if 0, then just plots summary of corelation coeff.
useCorr =1; % if 1, then cauclate pearson's corr, if 0, then euclid dist (5ms bins), to ask about similarity,

premotor_wind = [-0.03 0.02]; % for cross correlation
removeNoiseTrials = 0;

lt_neural_BatchSmth_CtxtSep(DATAllSwitches, MotifSets, premotor_wind, ...
    useCorr, plotRaw, removeNoiseTrials)


%% ==== for each channel, PLOT get deviation for each switch, distribution across switches

close all;

% ====== input params
motifstoplot = {'ab', 'jb', 'g'}; % if empty, plots all
motifstoplot = {'a(b)', 'j(b)', 'h(g)'}; % if empty, plots all
motifstoplot = {'a(b)', 'a(a)b', 'j(j)b', 'h(g)'}; % if empty, plots all
fs = 30000;


chanstoplot = find(~cellfun('isempty', DATAllSwitches.switch(1).motif(1).batchinorder(1).DatAll));
numswitches = length(DATAllSwitches.switch);
nummotifs = length(DATAllSwitches.switch(1).motif);
assert(length(DATAllSwitches.switch(1).motif(1).batchinorder)==2, 'note, this code wont work since tries to take diference');


for cc = chanstoplot
    
    figcount=1;
    subplotrows=4;
    subplotcols=6;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    for mm = 1:nummotifs
        
        motif = DATAllSwitches.switch(1).motif(mm).batchinorder(1).motifname;
        if ~isempty(motifstoplot)
            if ~any(ismember(motifstoplot, motif))
                continue
            end
        end
        
        
        DprimeAll = []; % one collect, one row for each switch
        
        for i=1:numswitches
            
            DatMeans ={};
            DatMedians = {};
            DatSEMs = {};
            CondAll = {};
            DatSTDs = {};
            for bb = 1:length(DATAllSwitches.switch(i).motif(mm).batchinorder)
                % i.e., each batch corresponds to a "condition" (e.g.
                % WNon/off; or DIR/UNDIR);
                datmat = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAll{cc};
                t = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).t;
                FF = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).FF;
                cond = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).condition;
                motif = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).motifname;
                
                try
                    pretime = DATAllSwitches.switch(i).params.pretime;
                    pretime = pretime/fs;
                catch err
                    disp('NOTE!!! making pretime 0.1');
                    pretime = 0.1;
                end
                
                % --- [INDIVIDUAL TRIALS] mean FR in a premotor window
                
                
                % --- [STATS FOR THIS SWITCH]
                datmedian = median(datmat,1);
                datmean = mean(datmat,1);
                datSTD = std(datmat);
                datSEM = lt_sem(datmat);
                
                
                % ============ save
                DatMeans = [DatMeans datmean];
                DatSEMs = [DatSEMs datSEM];
                CondAll = [CondAll cond];
                DatSTDs = [DatSTDs datSTD];
                DatMedians = [DatMedians datmedian];
            end
            
            % ----------- take difference of means, scaled by standard
            % deviation (d')
            numer = (DatMeans{2}-DatMeans{1});
            denom = sqrt(0.5*(DatSTDs{1}.^2 + DatSTDs{2}.^2));
            dprime = numer./denom;
            
            medianDiff = DatMedians{2} - DatMedians{1};
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['ch' num2str(cc) '-sw' num2str(i) '-' motif '[' CondAll{2} 'minus' CondAll{1} ']']);
            ylabel('dprime');
            
            plot(t, dprime, 'k');
            
            line([pretime pretime],ylim);
            lt_plot_zeroline;
            axis tight;
            
            % ================= collect
            DprimeAll = [DprimeAll; dprime];
            
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title('mean of abs');
        shadedErrorBar(t, mean(abs(DprimeAll),1), lt_sem(abs(DprimeAll)), {'Color', 'r'}, 1);
        
        % ======== mean for this motif, across all switches
        
        linkaxes(hsplots, 'xy');
    end
end



%% ====== INDIVIDUAL TRIALS
close all;
motifstoplot = {'a(b)', 'j(b)', '(a)ab', 'h(g)'}; % if empty, plots all
fs = 30000;
premotor_window = [-0.035 0.025]; % rel onset

chanstoplot = find(~cellfun('isempty', DATAllSwitches.switch(1).motif(1).batchinorder(1).DatAll));
numswitches = length(DATAllSwitches.switch);
nummotifs = length(DATAllSwitches.switch(1).motif);
assert(length(DATAllSwitches.switch(1).motif(1).batchinorder)==2, 'note, this code wont work since tries to take diference');


for cc = chanstoplot
    
    figcount=1;
    subplotrows=4;
    subplotcols=6;
    fignums_alreadyused=[];
    hfigs=[];
    for mm = 1:nummotifs
        
        motif = DATAllSwitches.switch(1).motif(mm).batchinorder(1).motifname;
        if ~isempty(motifstoplot)
            if ~any(ismember(motifstoplot, motif))
                continue
            end
        end
        
        
        FrateAllAll = []; % all trials for this motif and chan (i.e. across all switches)
        PitchAllAll = [];
        hsplots1 = [];
        hsplots2 = [];
        
        for i=1:numswitches
            
            CondAll = {}; % one entry for each block
            condmat = []; % one entry fro each trial
            FrateAllTrials = [];
            PitchAllTrials = [];
            
            for bb = 1:length(DATAllSwitches.switch(i).motif(mm).batchinorder)
                % i.e., each batch corresponds to a "condition" (e.g.
                % WNon/off; or DIR/UNDIR);
                datmat = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAll{cc};
                t = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).t;
                FF = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).FF;
                cond = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).condition;
                motif = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).motifname;
                
                try
                    pretime = DATAllSwitches.switch(i).params.pretime;
                    pretime = pretime/fs;
                catch err
                    disp('NOTE!!! making pretime 0.1');
                    pretime = 0.1;
                end
                
                % --- [INDIVIDUAL TRIALS] mean FR in a premotor window
                windwind = premotor_window+pretime;
                indtmp = t>=windwind(1) & t<windwind(2);
                FrateAllTrials = [FrateAllTrials; mean(datmat(:,indtmp),2)];
                PitchAllTrials = [PitchAllTrials; FF];
                condmat = [condmat; bb*ones(size(FF))];
                
                CondAll = [CondAll cond]; % length 2
                
                % -------------- collect across all switches
                FrateAllAll = [FrateAllAll; mean(datmat(:,indtmp),2)]; % all trials for this motif and chan (i.e. across all switches)
                PitchAllAll = [PitchAllAll; FF];
                
            end
            
            % =================== 1) PLOT regression, Frate vs. FF (including
            % both pre and post switch)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots1 = [hsplots1 hsplot];
            ylabel('Frate vs. Pitch');
            title(['ch' num2str(cc) '-sw' num2str(i) '-' motif '[' CondAll{2} '(r)-' CondAll{1} '(k)]']);
            hold on;
            indstmp = condmat==1;
            plotcol = 'k';
            plot(PitchAllTrials(indstmp), FrateAllTrials(indstmp), 'o', 'Color', plotcol);
            indstmp = condmat==2;
            plotcol = 'r';
            plot(PitchAllTrials(indstmp), FrateAllTrials(indstmp), 'o', 'Color', plotcol);
            hold on;
            lt_regress(FrateAllTrials, PitchAllTrials, 0, 0, 1, 1, 'b', 1)
            axis tight
            
            
            % =================== 2) cross correaltion Frate and FF
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots2 = [hsplots2 hsplot];
            title('xcov, left is frate leads');
            xlabel('frate <--> pitch');
            windowmax = 10;
            binsize=1;
            [xx, lags] = xcov(FrateAllTrials, PitchAllTrials, ceil(windowmax/binsize), 'coeff');
            plot(lags*binsize, xx, '-k');
            lt_plot_zeroline;
        end
        
        % ======== across all switches
        % ------------- 1) regression
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots1 = [hsplots1 hsplot];
        ylabel('Frate vs. Pitch');
        title('ALL SWITCHES');
        lt_regress(FrateAllAll, PitchAllAll, 1, 0, 1, 1, 'b', 1)
        axis tight
        
        % --------------- 2) xcov
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots2 = [hsplots2 hsplot];
        title('ALL TRIALS');
        xlabel('frate <--> pitch');
        windowmax = 10;
        binsize=1;
        [xx, lags] = xcov(FrateAllAll, PitchAllAll, ceil(windowmax/binsize), 'coeff');
        plot(lags*binsize, xx, '-k', 'LineWidth', 2);
        lt_plot_zeroline;
        
        linkaxes(hsplots1, 'xy');
        linkaxes(hsplots2, 'xy');
        
    end
    
end


%% ======================== CALCULATE CROSS COV BETWEEN BRAIN REGIONS
close all;
motifstoplot = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)', 'jbh(h)', 'g(h)'};
windowmax = 0.04;
binsize = 0.002;
premotorwind = [-0.045 0.025]; % use for ch14-21 on 11/12, morning.
fs = 30000;
plotRawDatOnly = 0;
chanstoplot = [14 21];
removenoise = 1;

lt_neural_BatchSmth_CrossRegion(DATAllSwitches, motifstoplot, chanstoplot, ...
    windowmax, binsize, premotorwind, fs, plotRawDatOnly, removenoise);



