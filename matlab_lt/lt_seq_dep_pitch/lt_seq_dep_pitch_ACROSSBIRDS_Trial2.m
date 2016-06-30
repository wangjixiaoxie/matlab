%% LT 6/16/16 - PLOT TRIAL-BY-TRIAL LEARNING, COMPARE TRAJECTORY RELATIVE TO LEARN MAGNITUDE, ETC (A LA SOBER)
function lt_seq_dep_pitch_ACROSSBIRDS_Trial2(SeqDepPitch_AcrossBirds, PARAMS, HourBins, plotCents)

% --- 1) replot a la sam sober (in cents)
% --- 2) what if plot to point when learning is at magnitude of shift from
% Sam's paper
% --- 3) Overnight vs. overday

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

OnlyExptsWithNoStartDelay=0; % removes expts with that start with delay from WN to learning.
DayWindow=[-3 6]; % [baseDays, WNdays];

% HourBins=[9:2:19]; % in 24 hour clock.
% HourBins=[8:3:20]; % in 24 hour clock.
minDataInBin=3; % ignore bin unless has at least this many points

if ~exist('plotCents', 'var');
    plotCents=0;
end

clear DATSTRUCT;

%% ==== Plot all birds like sam sober (but separate by syl, smooth across days, plot by rendition)

count=1;

for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ===== SKIP IF HAD DELAY FROM WN TO EXPT START
        if OnlyExptsWithNoStartDelay==1
           if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode>0
               disp(['SKIPPED (delay after WN start): ' birdname '-' exptname]);
               continue
           end
            
        end
        
        
        % ==== other infor about this expt
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        NumDays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SingleDirLastDay;
        FirstDay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
        
        for j=1:length(SylsUnique)
            
            syl=SylsUnique{j};
            % ===== extract data for this syl
            ffvals_alldays=[];
            tvals_alldays=[];
            
            for day=1:NumDays;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                    % -- then use time window
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
                    tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                else
                    % -- then use all day dat
                    ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day});
                    tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                end
                
                ffvals_alldays=[ffvals_alldays ffvals];
                tvals_alldays=[tvals_alldays tvals];
            end
            
            % -- targdir sign if needed
            targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            
            
            % ---- convert tvals to units of days rel to start
            tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, tvals_alldays);
            tmp.FinalValue;
            if floor(tmp.FinalValue(1))~=1;
                disp('PROBLEM, WHY IS FIRST DAY NOT DAY 1?');
            end
            
            % -- base mean
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                baseMeanFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
            else
                baseMeanFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            end
            
            DATSTRUCT.data(count).syl(j).syl=syl;
            DATSTRUCT.data(count).syl(j).baseMeanFF=baseMeanFF;
            DATSTRUCT.data(count).syl(j).ffvals_alldays=ffvals_alldays;
            DATSTRUCT.data(count).syl(j).tvals_alldays=tvals_alldays;
            DATSTRUCT.data(count).syl(j).dayVals_alldays=tmp.FinalValue;
            DATSTRUCT.data(count).syl(j).same=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            DATSTRUCT.data(count).syl(j).targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            
            DATSTRUCT.data(count).WNday1=WNday1;
            DATSTRUCT.data(count).SingleDirLastDay=NumDays;
            DATSTRUCT.data(count).birdname=birdname;
            DATSTRUCT.data(count).exptname=exptname;
            DATSTRUCT.data(count).targdir=targdir;
        end
        
        count=count+1;
        
    end
    
end

NumExpt=length(DATSTRUCT.data);


%% ======= PLOT MEAN LEARNING TRAJECTORY
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

plotXbinNum=1; % 1 uses bin num, 0 uses days;

Yall_same=[];
Yall_diff=[];
Yall_targ=[];

for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    
    % ==== PLOT FOR THIS EXPT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    
    
    numsyls=length(DATSTRUCT.data(i).syl);
    
    % ---- for means across expt
    alldat_x_same=[];
    allday_y_same=[];
    alldat_x_diff=[];
    allday_y_diff=[];
    
    for ii=1:numsyls
        
        ffvals=DATSTRUCT.data(i).syl(ii).ffvals_alldays;
        dayVals=DATSTRUCT.data(i).syl(ii).dayVals_alldays;
        same=DATSTRUCT.data(i).syl(ii).same;
        targ=DATSTRUCT.data(i).syl(ii).targ;
        WNday1=DATSTRUCT.data(i).WNday1;
        baseMeanFF=DATSTRUCT.data(i).syl(ii).baseMeanFF;
        
        if targ==1
            color='k';
        elseif same==1 & targ==0
            color='b';
        elseif same==0
            color='r';
        end
        
        % --- extract data for timepoints surrounding WN start
        minday=WNday1+DayWindow(1);
        maxday=WNday1+DayWindow(2)-1;
        inds=floor(dayVals)>=minday & floor(dayVals)<=maxday;
        
        ffvals=ffvals(inds)-baseMeanFF;
        dayVals=dayVals(inds);
        
        
        % - flip ffvals sign if needed
        ffvals=DATSTRUCT.data(i).targdir.*ffvals;
        
        % ========== CONVERT TO CENTS?
        if plotCents==1
            ffvals=1200*log2((baseMeanFF+ffvals)/baseMeanFF);
        end
        
        % ----- BIN BY SUB-DAY PORTIONS
        firstday=floor(min(dayVals)); % for the extracted data
        lastday=floor(max(dayVals));
        
        % -- construct all bin edges across all days
%         tmp=repmat([firstday:lastday]', 1, length(HourBins));
%         tmp2=repmat(HourBins, lastday-firstday+1, 1)./24;
        tmp=repmat([minday:maxday]', 1, length(HourBins));
        tmp2=repmat(HourBins, maxday-minday+1, 1)./24;
        
        DayValBinEdges=tmp+tmp2;
        
        % --- Go through each day - for each day go through all bins
        %         DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(numel(DayValBinEdges)-size(DayValBinEdges,1))=[];
        bincount=0;
        for j=1:size(DayValBinEdges,1); % days
            for jj=1:size(DayValBinEdges, 2)-1;
                binedges=[DayValBinEdges(j, jj) DayValBinEdges(j, jj+1)];
                
                bincount=bincount+1;
                binCenter=mean(binedges);
                % ===== COLLECT DATA
                inds=dayVals>binedges(1) & dayVals<=binedges(2);
                
                % ------ OUTPUT
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).binCenter=binCenter;
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).binEdges=binedges;
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).binNumIncludEmpty=bincount;
                
                if sum(inds)<minDataInBin
                    % skip any bins with less
                    
                    % --- save some info
                    DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).meanFF=nan;
                    DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).ffstd=nan;
                    DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).ffsem=nan;
                    DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).N=nan;
                    DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).dayvalsmean=nan;
                    
                    continue
                end
                
                
                fftemp=ffvals(inds);
                ffmean=mean(fftemp);
                ffstd=std(fftemp);
                ffsem=lt_sem(fftemp);
                
                dayvalstemp=dayVals(inds);
                dayvalsmean=mean(dayvalstemp);
                
                % ---- OUTPUT
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).meanFF=ffmean;
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).ffstd=ffstd;
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).ffsem=ffsem;
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).N=numel(fftemp);
                DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(bincount).dayvalsmean=dayvalsmean;
            end
            
        end
        
        % ====== PLOT THIS SYL
        xday=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.dayvalsmean]; % day vals of each bin (mean)
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.meanFF]; % ff mean in each bin
        ysem=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.ffsem];
        xbinnum=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binNumIncludEmpty];
        
        if DATSTRUCT.data(i).syl(ii).targ==1;
            plotcol='k';
        elseif  DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            plotcol='b';
        elseif  DATSTRUCT.data(i).syl(ii).same==0; % diff type
            plotcol='r';
        end
        
        %         shadedErrorBar(xday, y, ysem, {'Color', color}, 1);
        if plotXbinNum==1
            if DATSTRUCT.data(i).syl(ii).targ==1;
                plot(xbinnum, y, '-', 'Color', plotcol, 'Linewidth', 2);
            else
                plot(xbinnum, y, '-', 'Color', plotcol);
            end
        else
            plot(xday, y, '-', 'Color', plotcol);
        end
        
        % === collect for means
        if plotXbinNum==1
            if DATSTRUCT.data(i).syl(ii).targ==1;
                % - across all expt
                Yall_targ=[Yall_targ; y];
                
            elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
                % - for this expt
                alldat_x_same=[alldat_x_same; xbinnum];
                allday_y_same=[allday_y_same; y];
                % - across all expt
                Yall_same=[Yall_same; y];
            elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
                alldat_x_diff=[alldat_x_diff; xbinnum];
                allday_y_diff=[allday_y_diff; y];
                
                % - across all expt
                Yall_diff=[Yall_diff; y];
            end
        end
        
    end
    
    % ===== PLOT MEANS
    % -- same
    if size(alldat_x_same,1)>1
        xmean=mean(alldat_x_same, 1);
        ymean=nanmean(allday_y_same, 1);
        
        plot(xmean, ymean, '-', 'Color', 'b','LineWidth', 2);
    end
    % -- diff
    if size(alldat_x_diff,1)>1
        xmean=mean(alldat_x_diff, 1);
        ymean=nanmean(allday_y_diff, 1);
        
        plot(xmean, ymean, '-', 'Color', 'r','LineWidth', 2);
    end
    
    % === put lines for day transitions
    if plotXbinNum==1
        % plot transitions between days as function of bin num
        NumBinsInDay=size(DayValBinEdges, 2)-1;
        NumDays=DayWindow(2)-DayWindow(1);
        transitions=NumBinsInDay+0.5:NumBinsInDay:NumBinsInDay*(NumDays-1)+0.5;
        for k=1:length(transitions)
            line([transitions(k) transitions(k)], ylim, 'Color','m');
        end
        % - for WN onset
        wnon=NumBinsInDay*(-1*DayWindow(1))+0.5;
        line([wnon wnon], ylim, 'Color','c')
    else
        line([WNday1 WNday1], ylim, 'Color','m');
    end
    line(xlim, [0 0], 'Color','m');
    
    % -- line for ~80 hz learning
    line(xlim, [80 80], 'Color', 'k','LineStyle', '--');
end


% ===== PLOT MEAN FOR ALL EXPERIMENTS
lt_figure; hold on;
if plotXbinNum==1
% -- targ
Ymean=nanmean(Yall_targ);
Ysem=lt_sem(Yall_targ);

shadedErrorBar(xbinnum, Ymean, Ysem,{'Color','k'}, 1);    
% -- same
Ymean=nanmean(Yall_same);
Ysem=lt_sem(Yall_same);

shadedErrorBar(xbinnum, Ymean, Ysem, {'Color','b'}, 1);    

% -- diff
Ymean=nanmean(Yall_diff);
Ysem=lt_sem(Yall_diff);

shadedErrorBar(xbinnum, Ymean, Ysem, {'Color','r'}, 1);    


        % plot transitions between days as function of bin num
        NumBinsInDay=size(DayValBinEdges, 2)-1;
        NumDays=DayWindow(2)-DayWindow(1);
        transitions=NumBinsInDay+0.5:NumBinsInDay:NumBinsInDay*(NumDays-1)+0.5;
        for k=1:length(transitions)
            line([transitions(k) transitions(k)], ylim, 'Color','m');
        end
        % - for WN onset
        wnon=NumBinsInDay*(-1*DayWindow(1))+0.5;
        line([wnon wnon], ylim, 'Color','c')

    line(xlim, [0 0], 'Color','m');
    
    % -- line for ~80 hz learning
    line(xlim, [80 80], 'Color', 'k','LineStyle', '--');
end

%% ==== note
lt_figure;
lt_plot_annotation(1, ['NOTE: BINNED meanFF sign already flipped'], 'r')

%% ====== PLOT SCATTER OF TARG SHIFT VS. NONTARG SHIFT - DURING WN WHEN TARG IS AROUND 
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    hsplots=[];
    
    numsyls=length(DATSTRUCT.data(i).syl);
    WNday1=DATSTRUCT.data(i).WNday1;
    
    ytarg=[];
    ysame=[];
    ydiff=[];    
    
    % ======= COLLECT targ, same type, diff type (BASE VS. DURING LEARNING)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++= BASE    
    for ii=1:numsyls
        
        inds=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]<WNday1);
        
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(inds).meanFF]; % ff across all bins
        
        % - OUTPUT
        if DATSTRUCT.data(i).syl(ii).targ==1; % target
            ytarg=[ytarg; y];
        elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            ysame=[ysame; y];
        elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
            ydiff=[ydiff; y];
        end
    end
    
    % ==== PLOT FOR THIS EXPT
    % ----------------- same
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    title(['BASELINE (same): ' birdname '-' exptname]);
    xlabel('target');
    ylabel('nontarg');
    
    for j=1:size(ysame,1);
        plot(ytarg, ysame(j, :), 'bo');
        
    end
    % --- mean
    ymean=nanmean(ysame,1);
    if ~isempty(ymean)
        plot(ytarg, ymean, 'bo', 'MarkerFaceColor', 'b');
    end
    lt_plot_zeroline; lt_plot_zeroline_vert;
    
    % -------------------- diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
title(['BASE (diff): ' birdname '-' exptname]);
    xlabel('target');
    ylabel('nontarg');
    for j=1:size(ydiff,1);
        plot(ytarg, ydiff(j, :), 'ro');
    end
    % --- mean
    ymean=nanmean(ydiff,1);
    if ~isempty(ymean)
        plot(ytarg, ymean, 'ro', 'MarkerFaceColor', 'r');
    end
        lt_plot_zeroline; lt_plot_zeroline_vert;
% +++++++++++++++++++++++++++++++++++++++++++++++++++++

    ytarg=[];
    ysame=[];
    ydiff=[];    

    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++ WN ON    
    for ii=1:numsyls
        
        inds=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]>=WNday1);
        
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(inds).meanFF]; % ff across all bins
        
        % - OUTPUT
        if DATSTRUCT.data(i).syl(ii).targ==1; % target
            ytarg=[ytarg; y];
        elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            ysame=[ysame; y];
        elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
            ydiff=[ydiff; y];
        end
    end
    
    % ==== PLOT FOR THIS EXPT
    % -- same
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
title(['WN ON (same): ' birdname '-' exptname]);
    xlabel('target');
    ylabel('nontarg');
    
    for j=1:size(ysame,1);
        plot(ytarg, ysame(j, :), 'bo');
    end
    % --- mean
    ymean=nanmean(ysame,1);
    if ~isempty(ymean)
        plot(ytarg, ymean, 'bo', 'MarkerFaceColor', 'b');
    end
    lt_plot_zeroline; lt_plot_zeroline_vert;
    
    % -- diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
title(['WN ON (diff): ' birdname '-' exptname]);
    xlabel('target');
    ylabel('nontarg');
    for j=1:size(ydiff,1);
        plot(ytarg, ydiff(j, :), 'ro');
    end
    % --- mean
    ymean=nanmean(ydiff,1);
    if ~isempty(ymean)
        plot(ytarg, ymean, 'ro', 'MarkerFaceColor', 'r');
    end
        lt_plot_zeroline; lt_plot_zeroline_vert;
% +++++++++++++++++++++++++++++++++++++++++++++++++++++
    
linkaxes(hsplots, 'xy');
xlim([-200 200]), ylim([-200 200]);



end



%% === ACROSS ALL EXPERIMENTS
hsplots=[];
binsize=100;

lt_figure; hold on;
% ++++++++++++++++++++++++++++++++++++++++++ BASELINE [SAME]
hsplot=lt_subplot(2,2,1); hold on;
hsplots=[hsplots hsplot];
title(['BASELINE (SAME) (all expt)']);
xlabel('target');
ylabel('nontarg');

% --- TO COLLECT ACROSS ALL EXPT
Yall_targ=[];
Yall_same=[];
for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    
    numsyls=length(DATSTRUCT.data(i).syl);
    WNday1=DATSTRUCT.data(i).WNday1;
    
    ytarg=[];
    ysame=[];
    ydiff=[];
    
    % ======= COLLECT targ, same type, diff type (BASE VS. DURING LEARNING)
    for ii=1:numsyls
        
        inds=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]<WNday1);
        
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(inds).meanFF]; % ff across all bins
        
        % - OUTPUT
        if DATSTRUCT.data(i).syl(ii).targ==1; % target
            ytarg=[ytarg; y];
        elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            ysame=[ysame; y];
        elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
            ydiff=[ydiff; y];
        end
    end
    
    for j=1:size(ysame,1);
        plot(ytarg, ysame(j, :), 'b.');
        
        % --- COLLECT
        Yall_targ=[Yall_targ ytarg];
        Yall_same=[Yall_same ysame(j, :)];
    end
    
end

% --- plot smoothed
[~, inds]=sort(Yall_targ);
Yall_targ=Yall_targ(inds);
tmp=lt_running_stats(Yall_targ, binsize);
Yall_targ=tmp.Mean;
Yall_targ_sem=tmp.SEM;

Yall_same=Yall_same(inds);
tmp=lt_running_stats(Yall_same, binsize);
Yall_same=tmp.Mean;
Yall_same_sem=tmp.SEM;

shadedErrorBar(Yall_targ, Yall_same, Yall_same_sem, {'Color', 'k'}, 1);

% lt_regress(Yall_same, Yall_targ, 1, 0, 1, 1, 'b');
    lt_plot_zeroline; lt_plot_zeroline_vert;

    
% ++++++++++++++++++++++++++++++++++++++++++ BASELINE [DIFF]
hsplot=lt_subplot(2,2,2); hold on;
hsplots=[hsplots hsplot];
title(['BASELINE (DIFF) (all expt)']);
xlabel('target');
ylabel('nontarg');

% --- TO COLLECT ACROSS ALL EXPT
Yall_targ=[];
Yall_same=[];
Yall_diff=[];
for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    
    numsyls=length(DATSTRUCT.data(i).syl);
    WNday1=DATSTRUCT.data(i).WNday1;
    
    ytarg=[];
    ysame=[];
    ydiff=[];
    
    % ======= COLLECT targ, same type, diff type (BASE VS. DURING LEARNING)
    for ii=1:numsyls
        
        inds=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]<WNday1);
        
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(inds).meanFF]; % ff across all bins
        
        % - OUTPUT
        if DATSTRUCT.data(i).syl(ii).targ==1; % target
            ytarg=[ytarg; y];
        elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            ysame=[ysame; y];
        elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
            ydiff=[ydiff; y];
        end
    end
    
    for j=1:size(ydiff,1);
                plot(ytarg, ydiff(j, :), 'r.');

        % --- COLLECT
        Yall_targ=[Yall_targ ytarg];
        Yall_diff=[Yall_diff ydiff(j, :)];
    end
    
end

% --- plot smoothed
[~, inds]=sort(Yall_targ);
Yall_targ=Yall_targ(inds);
tmp=lt_running_stats(Yall_targ, binsize);
Yall_targ=tmp.Mean;
Yall_targ_sem=tmp.SEM;

Yall_diff=Yall_diff(inds);
tmp=lt_running_stats(Yall_diff, binsize);
Yall_diff=tmp.Mean;
Yall_diff_sem=tmp.SEM;

shadedErrorBar(Yall_targ, Yall_diff, Yall_diff_sem, {'Color', 'k'}, 1);

% lt_regress(Yall_diff, Yall_targ, 1, 0, 1, 1, 'r');
    lt_plot_zeroline; lt_plot_zeroline_vert;
    
    
    
    % ++++++++++++++++++++++++++++++++++++++++++ WN [SAME]
hsplot=lt_subplot(2,2,3); hold on;
hsplots=[hsplots hsplot];
title(['WN ON (SAME) (all expt)']);
xlabel('target');
ylabel('nontarg');

% --- TO COLLECT ACROSS ALL EXPT
Yall_targ=[];
Yall_same=[];
for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    
    numsyls=length(DATSTRUCT.data(i).syl);
    WNday1=DATSTRUCT.data(i).WNday1;
    
    ytarg=[];
    ysame=[];
    ydiff=[];
    
    % ======= COLLECT targ, same type, diff type (BASE VS. DURING LEARNING)
    for ii=1:numsyls
        
        inds=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]>=WNday1);
        
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(inds).meanFF]; % ff across all bins
        
        % - OUTPUT
        if DATSTRUCT.data(i).syl(ii).targ==1; % target
            ytarg=[ytarg; y];
        elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            ysame=[ysame; y];
        elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
            ydiff=[ydiff; y];
        end
    end
    
    for j=1:size(ysame,1);
        plot(ytarg, ysame(j, :), 'b.');
        
        % --- COLLECT
        Yall_targ=[Yall_targ ytarg];
        Yall_same=[Yall_same ysame(j, :)];
    end
    
end

% --- plot smoothed
[~, inds]=sort(Yall_targ);
Yall_targ=Yall_targ(inds);
tmp=lt_running_stats(Yall_targ, binsize);
Yall_targ=tmp.Mean;
Yall_targ_sem=tmp.SEM;

Yall_same=Yall_same(inds);
tmp=lt_running_stats(Yall_same, binsize);
Yall_same=tmp.Mean;
Yall_same_sem=tmp.SEM;

shadedErrorBar(Yall_targ, Yall_same, Yall_same_sem, {'Color', 'k'}, 1);

% lt_regress(Yall_same, Yall_targ, 1, 0, 1, 1, 'b');
    lt_plot_zeroline; lt_plot_zeroline_vert;
    
    
    % ++++++++++++++++++++++++++++++++++++++++++ WN ON [DIFF]
hsplot=lt_subplot(2,2,4); hold on;
hsplots=[hsplots hsplot];
title(['WN ON (DIFF) (all expt)']);
xlabel('target');
ylabel('nontarg');

% --- TO COLLECT ACROSS ALL EXPT
Yall_targ=[];
Yall_same=[];
Yall_diff=[];
for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    
    numsyls=length(DATSTRUCT.data(i).syl);
    WNday1=DATSTRUCT.data(i).WNday1;
    
    ytarg=[];
    ysame=[];
    ydiff=[];
    
    % ======= COLLECT targ, same type, diff type (BASE VS. DURING LEARNING)
    for ii=1:numsyls
        
        inds=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]>=WNday1);
        
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR(inds).meanFF]; % ff across all bins
        
        % - OUTPUT
        if DATSTRUCT.data(i).syl(ii).targ==1; % target
            ytarg=[ytarg; y];
        elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            ysame=[ysame; y];
        elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
            ydiff=[ydiff; y];
        end
    end
    
    for j=1:size(ydiff,1);
        plot(ytarg, ydiff(j, :), 'r.');
        
        % --- COLLECT
        Yall_targ=[Yall_targ ytarg];
        Yall_diff=[Yall_diff ydiff(j, :)];
    end
    
end

% --- plot smoothed
[~, inds]=sort(Yall_targ);
Yall_targ=Yall_targ(inds);
tmp=lt_running_stats(Yall_targ, binsize);
Yall_targ=tmp.Mean;
Yall_targ_sem=tmp.SEM;

Yall_diff=Yall_diff(inds);
tmp=lt_running_stats(Yall_diff, binsize);
Yall_diff=tmp.Mean;
Yall_diff_sem=tmp.SEM;

shadedErrorBar(Yall_targ, Yall_diff, Yall_diff_sem, {'Color', 'k'}, 1);


% lt_regress(Yall_diff, Yall_targ, 1, 0, 1, 1, 'r');
    lt_plot_zeroline; lt_plot_zeroline_vert;
    
    

% =======
linkaxes(hsplots, 'xy');


%% ======= LINEARLY STRETCH EACH EXPT - DOES EARLY LEARNING, MATCHED TO TARG LEARNING, SHOW NEGATIVE SHIFT OF DIFF TYPES?
% gets slope with intercept of one, least squares estimate

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

intercepts=[];
scaleFactor=[];
useBinNum=1; % otherwise will use day vals for regression
for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    numsyls=length(DATSTRUCT.data(i).syl);
    WNday1=DATSTRUCT.data(i).WNday1;
    
    % ==== FIGURE OUT HOW MUCH TO STRETCH/CONTRACT LEARNING
    ind=find([DATSTRUCT.data(i).syl.targ]==1); % targ syl
    firstWNbin=find([DATSTRUCT.data(i).syl(ind).BINNEDBYHOUR.binCenter]>=WNday1, 1,'first');
    
    ffvals=[DATSTRUCT.data(i).syl(ind).BINNEDBYHOUR(firstWNbin:end).meanFF];
    
    if useBinNum==1
        dayvals=[DATSTRUCT.data(i).syl(ind).BINNEDBYHOUR(firstWNbin:end).binNumIncludEmpty];
        dayvals=dayvals-firstWNbin+1; % first bin with dat is "1" (so preceding bin is 0);
    else
        dayvals=[DATSTRUCT.data(i).syl(ind).BINNEDBYHOUR(firstWNbin:end).dayvalsmean];
        dayvals=dayvals-DATSTRUCT.data(i).syl(ind).BINNEDBYHOUR(firstWNbin).binEdges(1); % intercept wanted to be ~0.
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
%     xlabel('day rel to first WN bin');
    ylabel('targ shift');
    
    
    % intercept constrained to 0 at the first
    indstmp=~isnan(ffvals);
    slope=ffvals(indstmp)/dayvals(indstmp);
    % slope=ffvals(indstmp)/[ones(1, sum(indstmp)); dayvals(indstmp)]; % lin
    % regression
    xtmp=[0 max(dayvals)];
    plot(dayvals, ffvals, 'ok');
    plot(xtmp, xtmp.*slope, '-k');
    
    lt_plot_annotation(1, ['slope=' num2str(slope)], 'r');
    
    if (0) % linear regres
        [~, bint, ~,~,~, sumstats]=lt_regress(ffvals, dayvals, 1, 0, 1, 1, 'k');
        % line([WNday1 WNday1], ylim, 'Color', 'm');
        intercepts=[intercepts sumstats.intercept];
    end
    
    DATSTRUCT.data(i).scaleFactor_Learn=slope;
    scaleFactor=[scaleFactor slope]; % = slope
    lt_plot_zeroline;
    ylim([0 200]);
end

lt_figure; hold on;
lt_plot_cdf(scaleFactor, 'b', 1);
xlabel('scale Factors (slopes)');
ylabel('frequency');
title('cdf');

minScaleFactor=input('minimum scale factor? ');


%% ======== WILL INTERPOLATE IN UPCOMING SECTION. WHAT X VALS TO USE FOR INTERPOLATE?

% ---- what xvals to use?
% METHOD 1 - use lowest maximum value (across all experiments). that way
% ensures gets all experiments, but throws out data after this value
% ============= WHAT IS MINIMUM(MAXIMUM(X))?
MinX_positive=min(Xall_targ(:,end));
MinX_negative=max(Xall_targ(:,1));
Xinterp=ceil(MinX_negative):floor(MinX_positive);

% METHOD 2 - ENTER BY HAND
Xinterp=-50:100;

%% ================================================= PLOT, SCALED BY SCALE [ALSO GET INTERPOLATED VALUES]
% FACTOR
% NOTE: bin 1 = first bin with WN. AFter got this, then scaled.
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

hsplots=[];
% plotXbinNum=1; % 1 uses bin num, 0 uses days;


Xall_targ=[];
Xall_same=[];
Xall_diff=[];
Yall_same=[];
Yall_diff=[];
Yall_targ=[];


ReplaceScaledValsWithInterpolated=1;

for i=1:NumExpt
    birdname=DATSTRUCT.data(i).birdname;
    exptname=DATSTRUCT.data(i).exptname;
    
    % --- SKIP IF SCALE FACTOR IS NEGATIVE OR VERY LOW
    if DATSTRUCT.data(i).scaleFactor_Learn<minScaleFactor
        continue
    end
    
    % ==== PLOT FOR THIS EXPT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    title([birdname '-' exptname]);
    
    numsyls=length(DATSTRUCT.data(i).syl);
    
    % ---- for means across expt
    alldat_x_same=[];
    allday_y_same=[];
    alldat_x_diff=[];
    allday_y_diff=[];
    
    for ii=1:numsyls
        
        same=DATSTRUCT.data(i).syl(ii).same;
        targ=DATSTRUCT.data(i).syl(ii).targ;
        WNday1=DATSTRUCT.data(i).WNday1;
        
        if targ==1
            color='k';
        elseif same==1 & targ==0
            color='b';
        elseif same==0
            color='r';
        end
                
        % ====== PLOT THIS SYL
        y=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.meanFF]; % ff mean in each bin
        ysem=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.ffsem];
        xbinnum=[DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binNumIncludEmpty];
        
        % - convert so that bin 1 is first WN bin
        firstWNbin=find([DATSTRUCT.data(i).syl(ii).BINNEDBYHOUR.binCenter]>=WNday1, 1,'first');
        xbinnum=xbinnum-firstWNbin+1;
        
        if DATSTRUCT.data(i).syl(ii).targ==1;
            plotcol='k';
        elseif  DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
            plotcol='b';
        elseif  DATSTRUCT.data(i).syl(ii).same==0; % diff type
            plotcol='r';
        end
        
        
        % ========== SCALE ALL BIN/TIME VALUES
        scaleFactor=DATSTRUCT.data(i).scaleFactor_Learn;
        xbinnum=xbinnum.*scaleFactor;
        
        
        % =================== INTERPOLATE
        if ReplaceScaledValsWithInterpolated==1
            % --- skip if has to extrapolate
            indstmp=~isnan(y);
            
            if max(xbinnum(indstmp)) < Xinterp(end) % largest x value with data
                disp(['CANT EXTRAPOLATE -' birdname '-' exptname]);
                continue;
            end
            Yinterp=interp1(xbinnum(indstmp), y(indstmp), Xinterp, 'spline', nan);
            
            xbinnum=Xinterp;
            y=Yinterp;
        end
        
            % ============== PLOT            
            if DATSTRUCT.data(i).syl(ii).targ==1;
                plot(xbinnum, y, '-', 'Color', plotcol, 'Linewidth', 2);
            else
                plot(xbinnum, y, '-', 'Color', plotcol);
            end
        
            
        % === collect for means
            if DATSTRUCT.data(i).syl(ii).targ==1;
                % - across all expt
                Xall_targ=[Xall_targ; xbinnum];
                Yall_targ=[Yall_targ; y];
                
            elseif DATSTRUCT.data(i).syl(ii).targ==0 &  DATSTRUCT.data(i).syl(ii).same==1; % same type
                % - for this expt
                alldat_x_same=[alldat_x_same; xbinnum];
                allday_y_same=[allday_y_same; y];
                % - across all expt
                Yall_same=[Yall_same; y];
                Xall_same=[Xall_same; xbinnum];
            elseif DATSTRUCT.data(i).syl(ii).same==0; % diff type
                alldat_x_diff=[alldat_x_diff; xbinnum];
                allday_y_diff=[allday_y_diff; y];
                
                % - across all expt
                Yall_diff=[Yall_diff; y];
                Xall_diff=[Xall_diff; xbinnum];
                
            end
        
    end
    
    % ===== PLOT MEANS
    % -- same
    if size(alldat_x_same,1)>1
        xmean=mean(alldat_x_same, 1);
        ymean=nanmean(allday_y_same, 1);
        
        plot(xmean, ymean, '-', 'Color', 'b','LineWidth', 2);
    end
    % -- diff
    if size(alldat_x_diff,1)>1
        xmean=mean(alldat_x_diff, 1);
        ymean=nanmean(allday_y_diff, 1);
        
        plot(xmean, ymean, '-', 'Color', 'r','LineWidth', 2);
    end
    
    % === put lines for day transitions
    % plot transitions between days as function of bin num
    NumBinsInDay=size(DayValBinEdges, 2)-1;
    NumDays=DayWindow(2)-DayWindow(1);
    transitions=NumBinsInDay+0.5:NumBinsInDay:NumBinsInDay*(NumDays-1)+0.5;
    transitions=(transitions - firstWNbin+1).*scaleFactor;
    for k=1:length(transitions)
        line([transitions(k) transitions(k)], ylim, 'Color','m');
    end
    % - for WN onset
    wnon=NumBinsInDay*(-1*DayWindow(1))+0.5;
    wnon=(wnon - firstWNbin+1)*scaleFactor;
    line([wnon wnon], ylim, 'Color','c')
    % --
    line(xlim, [0 0], 'Color','m');
    
    % -- line for ~80 hz learning
    line(xlim, [80 80], 'Color', 'k','LineStyle', '--');
end

linkaxes(hsplots, 'xy');




%% ========= PLOT [FOR BOTH INTERP AND NON INTERPOLATED DATA]
% ===== PLOT MEAN FOR ALL EXPERIMENTS
lt_figure; hold on;
lt_subplot(2,2,1); hold on;
title('all syls all expts, each plotted');
% --- targ
plot(Xall_targ', Yall_targ', '-k');
% --- same
plot(Xall_same', Yall_same', '-b');
% --- diff
plot(Xall_diff', Yall_diff', '-r');



%% ============================== IF INTERPOLATED, THEN CAN PLOT COLLAPSING ACROSS X VALUES
if ReplaceScaledValsWithInterpolated==1
    
    lt_subplot(2,2,2); hold on;
    title('mean, SEM [interpolated only]');
    
    % --- targ
    shadedErrorBar(Xall_targ(1,:), mean(Yall_targ, 1), lt_sem(Yall_targ), {'Color', 'k'}, 1);
     % --- same
    shadedErrorBar(Xall_same(1,:), mean(Yall_same, 1), lt_sem(Yall_same), {'Color', 'b'}, 1);
     % --- diff
    shadedErrorBar(Xall_diff(1,:), mean(Yall_diff, 1), lt_sem(Yall_diff), {'Color', 'r'}, 1);
   
end


%% ======== OTHER STUFF/PLOTS


% =============================== GET AVERAGE BY BINNING (COMBINING ALL
% POINTS, FOR GIVEN SYL, WITHIN A BIN, INTO ONE VAL (NO PSEUDOREPLICATION)

% --- what is minimun bin size? make this the size of the largest interval
% I have in my dataset. all other syls will combined all points within a
% bin into one point.
MaxDiff=[];

for j=1:size(Xall_diff, 1)
    xvals=Xall_diff(j, :);
    yvals=Yall_diff(j, :);
    inds=~isnan(yvals) & xvals>0;
    
    xvals=xvals(inds);
    yvals=yvals(inds);
    
    tmp=max(diff(yvals)); % largest difference between adjacents points (in units of bins, prescaling) for this syl
    
    MaxDiff=max([MaxDiff tmp]);
        
end
    

tmp=diff(Xall_targ, 1, 2)
diff(tmp, 1, 2)


% ============= WHAT IS MINIMUM(MAXIMUM(X))?
MinX=[];

for j=1:size(Xall_targ, 1)
    xvals=Xall_targ(j, :);

    
end




% TEMPORARY - CHOOSE A SCALED X VALUE, GET MEAN PITCH AT THAT TIME. 
BinsToAdd=2; %  =/- from the closest bin.
yall=[];
for j=1:size(Xall_diff, 1)
    xvals=Xall_diff(j, :);
    yvals=Yall_diff(j, :);
    inds=~isnan(yvals);
    xvals=xvals(inds);
    yvals=yvals(inds);
    
    
    [~, ind]=min(abs(xvals-65)); % timepoint of learning to check
    if length(yvals)>=ind+BinsToAdd
    ymean=yvals(ind-BinsToAdd:ind+BinsToAdd);
    else
    ymean=yvals(ind-BinsToAdd:ind);
    end        
    ymean=nanmean(ymean);
    
    yall=[yall ymean];
end





    
    
