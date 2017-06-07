function lt_seq_dep_pitch_ACROSSBIRDS_LearnVsNumCtxts(SeqDepPitch_AcrossBirds)
%% TO DO
% for learning rate, currently only keeping experiments with all 4 days
% having data. Modify it to calculate rate within each day by getting pitch
% at start and end of day, and dividing by numbher of rends. Doing this
% also allows me to keep even experiments that are missing some days.

%% lt 5/24/17 - number fo contexts greater, "drags" down learning in the target context?

NumBirds = length(SeqDepPitch_AcrossBirds.birds);
LearningOnDayNum = []; % if empty, then uses day 3-4 as normal. (i./e. day 1 is first WN day)
IgnorePu63ForLearningRate = 1; % keep this at 1, estimate for renditions is probably way off

%% === for each experiment, extract information to plot

NumSameTypes = [];
LearningTargHz = [];
LearningTargZscore = [];
LearningTargPercent = [];
TransProbTargContext = [];

LearningRateZ_all = [];
NumRends_all = [];

MaxNumBackSimilar_All = [];

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % ==== NUMBER OF SAME TYPES
        SametypeSyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
        %         disp(length(SametypeSyls));
        NumSameTypes = [NumSameTypes length(SametypeSyls)];
        
        % ==== OUT OF ALL SAME TYPES, WHAT IS MAXIMUM NUMBER-BACK?
        tmp = [];
        for j=1:length(SametypeSyls)
           syl = SametypeSyls{j};
           
           presylsimilar = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
           twosylsimilar = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ;
           
           contextsimilarity = 0;
           if presylsimilar==1 & twosylsimilar ==0
               contextsimilarity=1;
           elseif presylsimilar ==1 & twosylsimilar ==1
               contextsimilarity = 2;
           end
           
           tmp = [tmp; contextsimilarity];           
        end
        if ~isempty(tmp)
        MaxNumBackSimilar_All = [MaxNumBackSimilar_All max(tmp)];
        else
        MaxNumBackSimilar_All = [MaxNumBackSimilar_All nan];
        end
        
        % ==== LEARNING
        targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        targlearndir = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        if isempty(LearningOnDayNum)
            learninghz = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        else
            WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            numdaystoskip = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
            daytmp = WNday1+numdaystoskip+LearningOnDayNum-1;
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                learninghz = SeqDepPitch_AcrossBirds.birds{i}.experiment{...
                    ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(daytmp);
            else
                learninghz = SeqDepPitch_AcrossBirds.birds{i}.experiment{...
                    ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase(daytmp);
            end
            
            % -- aside, confirming that num skipped days is accurate
            assert(WNday1+numdaystoskip+3-1 == SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed(1), 'asdfasdf');
        end
        learninghz = learninghz*targlearndir;
        LearningTargHz = [LearningTargHz learninghz];
        
        
        % -- convert to z and percent
        basemeanFF = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).meanFF_RecalcDays;
        basestdFF = std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).rawFF_RecalcDays);
        learningZ = learninghz/basestdFF;
        learningPercent = learninghz/basemeanFF;
        
        LearningTargZscore = [LearningTargZscore learningZ];
        LearningTargPercent = [LearningTargPercent learningPercent];
        
        
        % --- FREQUENCY OF TARGET CONTEXT (RELATIVE TO OTHER CONTEXTS FOR
        % TARGET SYLLABLE)
        NumRendsOtherContexts = []; % from entire baseline
        for j=1:length(SametypeSyls)
            syl = SametypeSyls{j};
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                numrends = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
            else
                numrends = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF);
            end
            NumRendsOtherContexts = [NumRendsOtherContexts numrends];
        end
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
            NumRendsTargContext = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).rawFF_WithinTimeWindow);
        else
            NumRendsTargContext = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).rawFF);
        end
        
        transprob = NumRendsTargContext/(NumRendsTargContext + sum(NumRendsOtherContexts));
        TransProbTargContext = [TransProbTargContext transprob];
        
        %% FIGURE OUT LEARNING RATE
        NumRendUnlabeledSongs = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.LEARNING_RATE_AT_TARG.NumRendsUnlabeledSongs; % estimated number of renditions of the target
        if isempty(LearningOnDayNum) & sum(NumRendUnlabeledSongs>0)==4
            % then can do this, as taking day 3-4 of laerning
            learningRateZ = learningZ/sum(NumRendUnlabeledSongs);
            disp(sum(NumRendUnlabeledSongs));
            %            if sum(NumRendUnlabeledSongs)>10000
            %                keyboard
            %            end
            if IgnorePu63ForLearningRate==1
                if sum(NumRendUnlabeledSongs)>20000
                learningRateZ = nan;
                end
            end
            numrends = sum(NumRendUnlabeledSongs);
        else
            learningRateZ = nan;
            numrends = nan;
        end
        
        LearningRateZ_all = [LearningRateZ_all learningRateZ];
        NumRends_all = [NumRends_all numrends];
        
    end
end


%% === does max num back (out of all the nontargets) predict learning?

lt_figure; hold on;
lt_subplot(3,1,1); hold on;
xlabel('max contextual similarity');
lt_regress(LearningTargZscore, MaxNumBackSimilar_All, 1, 0, 1, 1, 'b',0);
xlim([-0.5 2.5])

%% ==== PLOT
lt_figure; hold on;

% =============== 1)
lt_subplot(3,1,1); hold on;
xlabel('num same type syls');
ylabel('learning at targ(hz)');

% plot(NumSameTypes, LearningTargHz, 'ok');
lt_regress(LearningTargHz, NumSameTypes, 1, 0, 1, 1, 'r', 0);
xlim([-0.5 4])

% =============== 2)
lt_subplot(3,1,2); hold on;
xlabel('num same type syls');
ylabel('learning at targ(zscore)');

% plot(NumSameTypes, LearningTargHz, 'ok');
lt_regress(LearningTargZscore, NumSameTypes, 1, 0, 1, 1, 'r', 0);
xlim([-0.5 4])

% =============== 2)
lt_subplot(3,1,3); hold on;
xlabel('num same type syls');
ylabel('learning at targ (percent)');

% plot(NumSameTypes, LearningTargHz, 'ok');
lt_regress(LearningTargPercent, NumSameTypes, 1, 0, 1, 1, 'r', 0);
xlim([-0.5 4])


% ============================= 2) USING TRANS PROB FOR TARG CONTEXT
lt_figure; hold on;

lt_subplot(4,1,1); hold on;
xlabel('num same type syl');
ylabel('targ context prob');
lt_regress(TransProbTargContext, NumSameTypes, 1, 0, 1, 1, 'r', 0);
xlim([-0.5 4])

lt_subplot(4,1,2); hold on;
xlabel('targ context prob');
ylabel('learning at targ (hz)');
lt_regress(LearningTargHz, TransProbTargContext, 1, 0, 1, 1, 'r', 0);
axis normal

lt_subplot(4,1,3); hold on;
xlabel('targ context prob');
ylabel('learning at targ (zscore)');
% inds = TransProbTargContext<1
lt_regress(LearningTargZscore, TransProbTargContext, 1, 0, 1, 1, 'r', 0);
axis normal


lt_subplot(4,1,4); hold on;
xlabel('targ context prob');
ylabel('learning at targ (%)');
lt_regress(LearningTargPercent, TransProbTargContext, 1, 0, 1, 1, 'r', 0);
axis normal



% ============================= 2) USING LAERNING RATE
lt_figure; hold on;

lt_subplot(5,1,1); hold on;
xlabel('num same type syl');
ylabel('num rends (estimated, across 4 days)');
lt_regress(NumRends_all, NumSameTypes, 1, 0, 1, 1, 'r', 0);
xlim([-0.5 4])

lt_subplot(5,1,2); hold on;
xlabel('trans prob targ context');
ylabel('num rends (estimated, across 4 days)');
lt_regress(NumRends_all, TransProbTargContext, 1, 0, 1, 1, 'r', 0);
axis normal

lt_subplot(5,1,3); hold on;
xlabel('num same type syl');
ylabel('learnig rate (z/rend)');
lt_regress(LearningRateZ_all, NumSameTypes, 1, 0, 1, 1, 'r', 0);
xlim([-0.5 4])

lt_subplot(5,1,4); hold on;
xlabel('targ context prob');
ylabel('learnig rate (z/rend)');
lt_regress(LearningRateZ_all, TransProbTargContext, 1, 0, 1, 1, 'r', 0);
axis normal

lt_subplot(5,1,5); hold on;
xlabel('targ context prob');
ylabel('learnig rate (z/rend)');
title('removing 100% targ context prob cases');
inds = TransProbTargContext<1;
lt_regress(LearningRateZ_all(inds), TransProbTargContext(inds), 1, 0, 1, 1, 'r', 0);
axis normal

% ============================== JUST COMPARE THOSE WITH OR WITHOUT OTHER
% CONTEXTS
lt_figure; hold on;

% ----------------- 1)
lt_subplot(2,1,1); hold on;
xlabel('1 context -- >1 context');
ylabel('learning Z')

Yall = {};

% - no contexts
inds = NumSameTypes==0;
X = 1;

Y = LearningTargZscore(inds);
plot(X,Y, 'ok');
lt_plot_bar(X, nanmean(Y), {'Errors', lt_sem(Y)});
Yall{X} = Y;


% - with other context
inds = NumSameTypes>0;
X = 2;

Y = LearningTargZscore(inds);
plot(X,Y, 'ok');
lt_plot_bar(X, nanmean(Y), {'Errors', lt_sem(Y)});
Yall{X} = Y;

% --
xlim([0 3]);

% -- rank sum
p = ranksum(Yall{1}, Yall{2});
lt_plot_pvalue(p, 'ranksum', 1);

% ----------------- 2)
lt_subplot(2,1,2); hold on;
xlabel('1 context -- >1 context');
ylabel('learning rate Z')

Yall = {};

% - no contexts
inds = NumSameTypes==0;
X = 1;

Y = LearningRateZ_all(inds);
plot(X,Y, 'ok');
lt_plot_bar(X, nanmean(Y), {'Errors', lt_sem(Y)});
Yall{X} = Y;


% - with other context
inds = NumSameTypes>0;
X = 2;

Y = LearningRateZ_all(inds);
plot(X,Y, 'ok');
lt_plot_bar(X, nanmean(Y), {'Errors', lt_sem(Y)});
Yall{X} = Y;

% --
xlim([0 3]);

% -- rank sum
p = ranksum(Yall{1}, Yall{2});
lt_plot_pvalue(p, 'ranksum', 1);
