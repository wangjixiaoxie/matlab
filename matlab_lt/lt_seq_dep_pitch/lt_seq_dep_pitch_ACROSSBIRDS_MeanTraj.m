function lt_seq_dep_pitch_ACROSSBIRDS_MeanTraj(SeqDepPitch_AcrossBirds, PARAMS, OnlyExptsWithNoStartDelay, TakeIntoAccountStartDelay, plotCents, DayWindow, throwOutIfAnyDayEmpty, plotZ, recalcBaseMean)
%% === notes
% plots relative to base extracted in these days.

%% PLOT MEAN TRAJECTORY OF LEARNING

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('plotCents', 'var');
    plotCents=0;
end



%% ==== Plot all birds like sam sober (but separate by syl, smooth across days, plot by rendition)

exptcount=1;

FFmeansTargALL=[];
ExptCountTargALL=[];

FFmeansSameALL=[];
ExptCountSameALL=[];

FFmeansDiffALL=[];
ExptCountDiffALL=[];

for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % ===== SKIP IF HAD DELAY FROM WN TO EXPT START
        if OnlyExptsWithNoStartDelay==1
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode>0
                disp(['SKIPPED (delay after WN start): ' birdname '-' exptname]);
                continue
            end
        end
        
        % ========== FIGURE OUT WHAT DAYS TO GET
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        
        % -- baseline days
        base1=WNday1+DayWindow(1);
        baseEnd=WNday1-1;
        % -- WN days
        if TakeIntoAccountStartDelay==1;
            daysDelay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
            if daysDelay>0
                disp(['ADDED ' num2str(daysDelay) ' days for: ' birdname '-' exptname]);
            end
        else
            daysDelay=0;
        end
        
        WN1=WNday1+daysDelay;
        WNend=WN1+DayWindow(2)-1;
        
        DaysToExtract=[base1:baseEnd WN1:WNend];
        BaseDays=[base1:baseEnd];
        WNDays=[WN1:WNend];
        % ================================================
        
        % ============= SKIP EXPT IF NOT ENOUGH DAYS
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        maxday=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase);
        if DaysToExtract(1)<1 | DaysToExtract(end)>maxday
            disp(['NOT ENOUGH DAYS, SKIPPING EXPT: ' birdname '-' exptname]);
            continue
        end
        
        % -- targdir sign if needed
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;

        % ----- EXTRACT FOR EACH SYL
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
%                 ffmeanvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(DaysToExtract);
                FFbase=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(BaseDays);
                FFwn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(WNDays);
                
                % - convert to actual ff
                tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                FFbase=FFbase+tmp;
                FFwn=FFwn+tmp;
                
                 baseSTD=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                 
            else
%                 ffmeanvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase(DaysToExtract)';
                FFbase=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF(BaseDays)';
                FFwn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF(WNDays)';
                
                baseSTD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
                basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            end
            
            % ============= get base mean, and convert to deviation from
            % mean
            if recalcBaseMean==1
                basemean=nanmean(FFbase);
            end
            
            ffmeanvals=[FFbase FFwn]-basemean;
            
            
            % ======= convert to cents?
            if plotZ==1
                
                ffmeanvals=(ffmeanvals)./baseSTD;
            else
                if plotCents==1
                    ffmeanvals=1200*log2((basemean+ffmeanvals)/basemean);
                end
            end
            
            % --- flip sign if needed
            ffmeanvals=ffmeanvals.*targdir;
            
            % ================= THROW OUT BECAUSE HAS AT LEAST A DAY WITH
            % NAN?
            if throwOutIfAnyDayEmpty==1
                if any(isnan(ffmeanvals))
                   disp(['THREW OUT SYL, has nan on some days: ' birdname '-' exptname '-' syl]);
                  continue 
               end
                
            end
            
            % ===== OUTPUT
            targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            same=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            
            if targ==1
                FFmeansTargALL=[FFmeansTargALL; ffmeanvals];
                ExptCountTargALL=[ExptCountTargALL exptcount];
                
            elseif targ==0 & same==1
                FFmeansSameALL=[FFmeansSameALL; ffmeanvals];
                ExptCountSameALL=[ExptCountSameALL exptcount];
                
            elseif same==0
                FFmeansDiffALL=[FFmeansDiffALL; ffmeanvals];
                ExptCountDiffALL=[ExptCountDiffALL exptcount];
                
            end
        end
        
        exptcount=exptcount+1;
    end
end
            
%% ============ PLOT GLOBAL MEAN

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

% =========== 1) MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean trajectory');
xlabel('day');
if plotCents==1
ylabel('shift (cents');
else
    ylabel('shift (hz)');
end

% - targ
plotcol='k';
y=nanmean(FFmeansTargALL,1);
ysem=lt_sem(FFmeansTargALL);
x=1:length(y);

shadedErrorBar(x, y, ysem, {'Color',plotcol}, 1);
% for each day, test if significantly diff from 0
for i=1:size(FFmeansTargALL, 2)
    ytmp=FFmeansTargALL(:, i);
    
    p=signrank(ytmp);
    if p<0.15;
        lt_plot_text(i, 1.1*nanmean(ytmp), ['p=' num2str(p, '%3.2g')], plotcol);
    end
end
% plot sample size
N=size(FFmeansTargALL,1);
lt_plot_annotation(1, ['N= ' num2str(N)], plotcol);


% - same
plotcol='b';
y=nanmean(FFmeansSameALL,1);
ysem=lt_sem(FFmeansSameALL);
x=1:length(y);

shadedErrorBar(x, y, ysem, {'Color',plotcol}, 1);
% for each day, test if significantly diff from 0
for i=1:size(FFmeansSameALL, 2)
    ytmp=FFmeansSameALL(:, i);
    
    p=signrank(ytmp);
    if p<0.15;
        lt_plot_text(i, 1.1*nanmean(ytmp), ['p=' num2str(p, '%3.2g')], plotcol);
    end
end
% plot sample size
N=size(FFmeansSameALL,1);
lt_plot_annotation(2, ['N= ' num2str(N)], plotcol);

% - diff
plotcol='r';
y=nanmean(FFmeansDiffALL,1);
ysem=lt_sem(FFmeansDiffALL);
x=1:length(y);

shadedErrorBar(x, y, ysem, {'Color',plotcol}, 1);
% for each day, test if significantly diff from 0
for i=1:size(FFmeansDiffALL, 2)
    ytmp=FFmeansDiffALL(:, i);
    
    p=signrank(ytmp);
    if p<0.15;
        lt_plot_text(i, 1.1*nanmean(ytmp), ['p=' num2str(p, '%3.2g')], plotcol);
    end
end

lt_plot_zeroline;
line([-DayWindow(1)+0.5, -DayWindow(1)+0.5], ylim)
% plot sample size
N=size(FFmeansDiffALL,1);
lt_plot_annotation(3, ['N= ' num2str(N)], plotcol);



% =========== 2) ALL TRAJS
lt_figure; hold on;
title('all trajectory');
xlabel('day');
if plotCents==1
ylabel('shift (cents');
else
    ylabel('shift (hz)');
end

% - targ
plotcol='k';
y=nanmean(FFmeansTargALL,1);
x=1:length(y);

plot(x, FFmeansTargALL, '-', 'Color', plotcol);


% - same
plotcol='b';
y=nanmean(FFmeansSameALL,1);
x=1:length(y);

plot(x, FFmeansSameALL, '-', 'Color', plotcol);


% - diff
plotcol='r';
y=nanmean(FFmeansDiffALL,1);
x=1:length(y);

plot(x, FFmeansDiffALL, '-', 'Color', plotcol);

%% ============= PLOT NONTARG SHIFT VS. TARG SHIFT
lt_figure; hold on;

NumExpts=max(ExptCountTargALL);

FFendSame=nanmean(FFmeansSameALL(:, [end-1 end]),2);
FFendTarg=nanmean(FFmeansTargALL(:, [end-1 end]),2);
FFendDiff=nanmean(FFmeansDiffALL(:, [end-1 end]),2);


% --------------------- SAME TYPE
lt_subplot(2,2,1); hold on;
title('same type (last 2 days)');
xlabel('targ shift (cents)');
ylabel('nontarg shift (cents)');
exptcountvals=ExptCountSameALL;
ffnontarg=FFendSame;
color='b';

X=[]; % targ
Y=[]; % nontarg;
for i=1:NumExpts
    
    % nontarg
    inds=exptcountvals==i;    
    yvals=ffnontarg(inds);
    
    % targ
    ind=ExptCountTargALL==i;
    xval=FFendTarg(ind);
    xvals=ones(1, length(yvals))*xval;
    
    % --- output
    X=[X xvals];
    Y=[Y yvals'];    
end

lt_regress(Y, X, 1, 0, 1, 1, color);    
lt_plot_zeroline; lt_plot_zeroline_vert;


% --------------------- DIFF TYPE
lt_subplot(2,2,2); hold on;
title('diff type');
xlabel('targ shift (cents)');
ylabel('nontarg shift (cents)');
exptcountvals=ExptCountDiffALL;
ffnontarg=FFendDiff;
color='r';

X=[]; % targ
Y=[]; % nontarg;
for i=1:NumExpts
    
    % nontarg
    inds=exptcountvals==i;    
    yvals=ffnontarg(inds);
    
    % targ
    ind=ExptCountTargALL==i;
    xval=FFendTarg(ind);
    xvals=ones(1, length(yvals))*xval;
    
    % --- output
    X=[X xvals];
    Y=[Y yvals'];    
end

lt_regress(Y, X, 1, 0, 1, 1, color);    
lt_plot_zeroline; lt_plot_zeroline_vert;








