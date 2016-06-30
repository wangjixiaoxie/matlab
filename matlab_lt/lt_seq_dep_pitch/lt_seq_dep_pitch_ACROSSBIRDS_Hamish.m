%% LT 6/20/16
function lt_seq_dep_pitch_ACROSSBIRDS_Hamish(SeqDepPitch_AcrossBirds, Params, plotExptRawDat)

%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

%%
useLearningRelLastBlineDay=0; % otherwise will use rel entire baseline
useWNday2=0; % otherwise will use 1. using 2 allows to get some expt without songs on day 1
takeDay2forExptWithNoDay1=1; % otherwise will throw out those experiments.



%% For each experiment, extract baseline LMAN bias, training direction, and metrics of learning
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


AFPbias_PBSminusMUSC_ALL=[];
TargDir_ALL=[];
Learning_day1_ALL=[];
CV_lastBaseDay_ALL=[];
CV_firstWNDay_ALL=[];
BirdNum_ALL=[];

for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        % ---- LMAN BIAS DURING BASELINE (using all baseline days)
        ffMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(targsyl).meanFF_WithinTimeWindow;
        ffPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).meanFF_WithinTimeWindow;
        
        AFPbias_PBSminusMUSC=ffPBS-ffMUSC;
        
        % ---- learning direction
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        
        % ---- learning on first day
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        WNday1actual=WNday1;
        if useWNday2==1
                       WNday1=WNday1+1;
        end
        
        % ===== FOR CERTAIN EXPERIMENTS, WITH no songs on day 1 (and
        % verified that did not sing, or sang very little, take day 2 as
        % start of learning
        if takeDay2forExptWithNoDay1==1
            if strcmp(birdname, 'pu53wh88') & strcmp(exptname, 'SeqDepPitchLMAN');
                disp('DID!!');
                WNday1=WNday1+1;
            end
            if strcmp(birdname, 'gr41gr90') & strcmp(exptname, 'SeqDepPitchLMAN');
                disp('DID!!');
                WNday1=WNday1+1;
            end
        end
        
        
        
        learning_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(WNday1);
        
        % ---- LEARNING as diff from first day to last bline day
        ffLastBlineDay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow{WNday1actual-1};
        ffFirstWNDay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow{WNday1};
        learning_day1_relLastBlineDay=mean(ffFirstWNDay)-mean(ffLastBlineDay);
        if useLearningRelLastBlineDay==1
            learning_day1=learning_day1_relLastBlineDay;
        end
        
        % ----- rd23 LM, use hand entered value using notcatch songs, since
        % not enough catch songs.
            if strcmp(birdname, 'rd23gr89') & strcmp(exptname, 'SeqDepPitchLMAN');
                % --- only do this if using learning on day 1
                if useWNday2==0
                learning_day1=3601.9-3569.3;
                end
            end
        
        
        % ---- CV on last baseline day
        ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow{WNday1actual-1};
        tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).Tvals_WithinTimeWindow{WNday1actual-1};
        %  plot, to visualize CV
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('blue = base; red = WN');
        plot(tvals, ffvals, 'ob');
        
        %        CV_lastBaseDay=std(ffvals)/mean(ffvals); OLD WAY, NOT USING CV.
        if ~isempty(ffvals);
            [~, ~, res]=lt_regress(ffvals, tvals, 0, 0, 0, 0, '');
            CV_lastBaseDay=std(res)./mean(ffvals);
        else
            CV_lastBaseDay=nan;
        end
        
        % ---- CV on first learning day
        ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow{WNday1};
        tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).Tvals_WithinTimeWindow{WNday1};
        plot(tvals-1, ffvals, 'or');
        
        %        CV_firstWNDay=std(ffvals)/mean(ffvals);
        if ~isempty(ffvals);
            [~, ~, res]=lt_regress(ffvals, tvals, 0, 0, 0, 0, '');
            CV_firstWNDay=std(res)./mean(ffvals);
        else
            CV_firstWNDay=nan;
        end
        
        % if CV_firstWNDay>CV_lastBaseDay
        %     disp(CV_firstWNDay);
        %     disp(CV_lastBaseDay);
        %     keyboard
        % end
        
        % ============== COLLECT
        disp(' ');
        disp(['---- ' birdname '-' exptname '-' targsyl]);
        disp(['targdir=' num2str(targdir)]);
        disp(['AFPbias=' num2str(AFPbias_PBSminusMUSC)]);
        disp(['learning_day1=' num2str(learning_day1)]);
        
        AFPbias_PBSminusMUSC_ALL=[AFPbias_PBSminusMUSC_ALL AFPbias_PBSminusMUSC];
        TargDir_ALL=[TargDir_ALL targdir];
        Learning_day1_ALL=[Learning_day1_ALL learning_day1];
        CV_lastBaseDay_ALL=[CV_lastBaseDay_ALL CV_lastBaseDay];
        CV_firstWNDay_ALL=[CV_firstWNDay_ALL CV_firstWNDay];
        BirdNum_ALL=[BirdNum_ALL i];
        
        % ====== PLOT RAW DAT FOR THIS DAY TO COMPARE TO EXTRACTED STATS
        if plotExptRawDat==1
            plotLarge=1;
            BirdToPlot=birdname;
            ExptToPlot=exptname;
            SylsToPlot={targsyl};
            overlayMeans=1;
            UseSylColors=0; % 0 is black;
            flipsign=1; % plot pos, if neg
            use_std=0; % applies to mean plots only!! (std, instead of sem)
            plotRawFF=1; % 1=raw FF, 0= baseline subtracted (MUSC-MUSC)
            OverlayLMANStats=1; % plots mean shift in single target learning window (defined below)
            OverlayMUSC_days=[];
            lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds, Params, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, OverlayMUSC_days, plotLarge)
            
            pause;
            close all;
        end
        
    end
end

%% ==== plot
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% --- WN up experiments
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('UP training');
inds=TargDir_ALL==1;

xlabel('afp bias (pbs - musc)');
ylabel('pitch shift');
afpbias=AFPbias_PBSminusMUSC_ALL(inds);
learning=Learning_day1_ALL(inds);
plot(afpbias, learning, 'ok');
lt_plot_zeroline; lt_plot_zeroline_vert;
xlim([-100 100]); ylim([-140 160]);

% --- WN down experiments
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DN training');
inds=TargDir_ALL==-1;

xlabel('afp bias (pbs - musc)');
ylabel('pitch shift');
afpbias=AFPbias_PBSminusMUSC_ALL(inds);
learning=Learning_day1_ALL(inds);
plot(afpbias, learning, 'ok');
lt_plot_zeroline; lt_plot_zeroline_vert;
xlim([-100 100]); ylim([-140 160]);


% ---- combined,
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all expt');
xlabel('AFP bias, relative to direction of training');
ylabel('pitch shift, in adaptive direction');

% plot all
    x=AFPbias_PBSminusMUSC_ALL.*TargDir_ALL;
    y=Learning_day1_ALL.*TargDir_ALL;
    lt_regress(y, x, 1, 0, 1, 1, 'k');

    % each bird diff col
    numbirds=max(unique(BirdNum_ALL));
plotcols=lt_make_plot_colors(numbirds, 0, 0);
for i=1:numbirds
    inds=BirdNum_ALL==i;
    x=AFPbias_PBSminusMUSC_ALL(inds).*TargDir_ALL(inds);
    y=Learning_day1_ALL(inds).*TargDir_ALL(inds);
    lt_plot(x, y, {'Color', plotcols{i}});
end
lt_plot_zeroline; lt_plot_zeroline_vert;
xlim([-100 100]); ylim([-140 160]);

% x=AFPbias_PBSminusMUSC_ALL.*TargDir_ALL;
% y=Learning_day1_ALL.*TargDir_ALL;
% lt_regress(y, x, 1, 0, 1, 1, 'k');
% lt_plot_zeroline; lt_plot_zeroline_vert;
% xlim([-100 100]); ylim([-140 160]);
% 


% ====================== DOES CV DECREASE FROM LAST BLINE DAY TO FIRST
% TRAINING DAY?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('pitcfh CV');
xlabel('last base day');
ylabel('first WN day');

lt_plot_45degScatter(CV_lastBaseDay_ALL, CV_firstWNDay_ALL, 'k', 1);
p=signrank(CV_lastBaseDay_ALL,CV_firstWNDay_ALL);
lt_plot_pvalue(p, 'signrank', 1);

% ========= DOES GREATER SHIFT IN DIRECTION OPPOSITE AFP BIAS CORRELATED
% WITH GREATER REDUCTION IN CV?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('pitch shift, relative direction of AFP bias');
ylabel('change in CV');

x=Learning_day1_ALL.*sign(AFPbias_PBSminusMUSC_ALL);
y=CV_firstWNDay_ALL-CV_lastBaseDay_ALL;
lt_regress(y, x, 1, 0, 1, 1, 'k');


% =========== IF SHIFT IS OPPOSITE AFP BIAS AND IN TARG DIR, SEE CV DOWN?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('cv reduced if goal is to shift in opposite dir?');
xlabel('afp bias, relative to training dir');
ylabel('CV change (diff)');

x=TargDir_ALL.*AFPbias_PBSminusMUSC_ALL; % 
y=CV_firstWNDay_ALL-CV_lastBaseDay_ALL;
lt_regress(y, x, 1, 0, 1, 1, 'k');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('cases where AFP bias was in same dir as training dir');
xlabel('pitch shift (adaptive dir)');
ylabel('CV change (diff)');

inds=sign(TargDir_ALL)==sign(AFPbias_PBSminusMUSC_ALL);
x=Learning_day1_ALL(inds).*TargDir_ALL(inds);
y=CV_firstWNDay_ALL(inds)-CV_lastBaseDay_ALL(inds);
lt_regress(y, x, 1, 0, 1, 1, 'b');
lt_plot_zeroline; lt_plot_zeroline_vert

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('cases where AFP bias was in opposite dir as training dir');
xlabel('pitch shift (adaptive dir)');
ylabel('CV change (diff)');

inds=sign(TargDir_ALL)~=sign(AFPbias_PBSminusMUSC_ALL);
x=Learning_day1_ALL(inds).*TargDir_ALL(inds);
y=CV_firstWNDay_ALL(inds)-CV_lastBaseDay_ALL(inds);
lt_regress(y, x, 1, 0, 1, 1, 'b');
lt_plot_zeroline; lt_plot_zeroline_vert











