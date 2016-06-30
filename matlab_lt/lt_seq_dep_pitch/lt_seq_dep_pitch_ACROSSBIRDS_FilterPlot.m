%% LT 4/29/15 - called by lt_seq_dep_pitch_ACROSSBIRDS
function [FILTERED_DATA, SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_FilterPlot(Params, SeqDepPitch_AcrossBirds)
% example inputs:
% Params.FilterPlot.sylID_filter{1}='presyl_similar_to_targ_presyl';
% Params.FilterPlot.sylID_filter{2}={0,1}; %
% 
% Params.FilterPlot.extra_filter=1; % if 1, then filter similar vs diff
% Params.FilterPlot.similar_only=1; % if 1, takes only similar, if 0, takes only different; only applies if extra_filter is on. 
% 
% % What days to plot for across days analysis?
% Params.FilterPlot.FirstDayToPlot=1 ; % lock to what day? (1st WN day would be 1);
% Params.FilterPlot.LastDayToPlot=[]; % lock end to what day? (5 would be 5th WN day, 'end' would be end of WN, [] would be openended);
% Params.FilterPlot.pre_days_to_plot=3; % dont change this.

% TO DO: N is incorrect for histogram


%% DEAULTS - other things - WOULD NOT NORMALLY CHANGE
skip_target=1; % if 1, skips collecting data for target. if 0 then collects
skip_exceptions=1; % if 1, skips exceptions entered by hand.

% Plot what data? (in progress)
data_to_plot='WN_start_daybins';
data_to_plot='WN_end_daybins';

% convert params to single variable names
sylID_filter=Params.FilterPlot.sylID_filter;
extra_filter=Params.FilterPlot.extra_filter;
similar_only=Params.FilterPlot.similar_only;

NumBirds_OneTarg=length(SeqDepPitch_AcrossBirds.birds);

Define_learning_as_consol_start=1; % other wise is end of WN





%% Filter data - i.e. pull out each datapoint 

FILTERED_DATA=[];
FILTERED_DATA.filter=sylID_filter{1};
for i=1:NumBirds_OneTarg;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:NumExperiments;
        syl_list=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % -- GET INFORMATION about the target syl
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
        if length(targsyl)>1;
            tmp=regexp(targsyl,'[A-Z]'); % find caps
            if tmp==1;
                targsyl_pre=nan;
            else
            targsyl_pre=targsyl(tmp-1);
            end
        end
        
        % learning metric for target
        LearnMetric_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        
        for iii=1:length(syl_list);
            syl=syl_list{iii};
            
            % -- PREPROCESSING ----------------------------------------
            
            % -- Skip certain ad hoc exceptions ------------------------
            if skip_exceptions==1;
                if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu11wh87') && strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID,'SyntaxDepPitchShift_cbDOWN');
                    if strcmp(syl,'dccB'); % skip because is also targeted by WN, but not designated as targsyl.
                        continue
                    end
                end
                
                if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu11wh87') && strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
                    if strcmp(syl,'bccB'); % skip because is also targeted by WN, but not designated as targsyl.
                        continue
                    end
                end
            end
            
            if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu37wh20') && strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
                if strcmp(syl,'bccB'); % skip because is also targeted by WN, but not designated as targsyl.
                    continue
                end
            end
            
            
            if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu64bk13');
                if strcmp(syl,'jB') || strcmp(syl,'jbB') || strcmp(syl,'jbbB') ; % skip because some get 100%Wn, some non-contingent WN (could not avoid).
                    continue
                end
            end
      
            % -----------------------------------------------
            
            % --- WHAT class is this syl?
            sylID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).(sylID_filter{1});
            
            % --- WHAT group is this data? ------------------------------
            % determining group is slightly different for different
            % dimensions
            group_num=nan;
            for j=1:length(sylID_filter{2}); % use for loop because not sure if is string or numerical
                if sylID_filter{2}{j}==sylID;
                    
                    % then is in group j
                    group_num=j;
                    break
                end
            end
            
            % warn user if did not find group
            if isnan(group_num);
                disp(['Warning -did not find group for bird ' num2str(i) '; experiment: ' num2str(ii) '; SKIPPING']);
                continue
            end
            
            
            % --- COLLECT data and group with others in this groupnum ---
            % baseline meanFF
            Baseline_FF_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            
            % 1) 3-day average to start and end WN
            WN_Start_Ind_Binned=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).FirstWNInd; % start of WN
            WN_start_daybins=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).rawFF{WN_Start_Ind_Binned}; % ff vals, last WN bin
            
            % 2) 3-day average to end WN
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1; % if there was only one targ, then last WN ind is straightforward
                WN_End_Ind_Binned=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).LastWNInd;
                WNend_daybins=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).rawFF{WN_End_Ind_Binned}; % ff vals,  1st WN bin
                
            elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0; % if started with one targ, then added another, then want last days of 1-targ, not the end of WN.
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot'); % inds should be here
                    snapshotfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
                    
                    WNend_daybins=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(snapshotfield{1}).(syl).FFvals;
                else
                    disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have - skipping'])
                    continue
                end
            end
            
                        
            % 3) Get day means for all days
            FF_minusbase_overdays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase;
            % - get indices for WN on and off
            WN_on_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            % - wn off day, depends on if one or one-->2 targ expt
            WN_off_day=nan;
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1; % if there was only one targ, then last WN ind is straightforward
                WN_off_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0; % if started with one targ, then added another, then want last days of 1-targ, not the end of WN.
                WN_off_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.DaysToMarkInds{1}; % assume this is the last one-targ WN day
            end
            if isnan(WN_off_day);
                disp('Problem WN off day not defined - skipping this experiment/bird');
                continue
            end

            
            % 5) GET LEARNING (relative to target)
            FFvals=WNend_daybins; % vals at end of WN
            FF_mean_minusbase=mean(FFvals)-Baseline_FF_mean; % converted to mean FF minus baseline
            
            % -- divide that value by pitch change at target.
            FF_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Learning_by_target.WN_end_daybins;
            Learning_rel_targ=FF_mean_minusbase/FF_targ; % learning as fraction of learning by target (minus respective baselines)

            
            % 6) Get consolidation data (start, end, change);
            CONSOLIDATION.FFminusbase_start = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_meanFF;
            CONSOLIDATION.FFminusbase_end = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolEnd_meanFF;
            
            CONSOLIDATION.FFsem_start = lt_sem(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_FFvals_minbase);
            CONSOLIDATION.FFsem_end = lt_sem(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolEnd_FFvals_minbase);
            
            
            Targ_start=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(targsyl).ConsolStart_meanFF;
            Targ_end=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(targsyl).ConsolEnd_meanFF;
            
            CONSOLIDATION.Learning_start = CONSOLIDATION.FFminusbase_start/Targ_start;
            CONSOLIDATION.Learning_end = CONSOLIDATION.FFminusbase_end/Targ_end;
            
            
            % 7) Define learning as consolidation start
            if Define_learning_as_consol_start==1;
            Learning_rel_targ=CONSOLIDATION.Learning_start;
            end
            
            % 8) Define learning using learning metric
            LearnMetric=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            LearnMetric_RelTarg=LearnMetric/LearnMetric_targ;
            
            
        % === PUT ALL DATA INTO OUTPUT STRUCTURE ---------------
            
            % if this is the target for this expt, put into target slot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                % then slide into target slot
                if isfield(FILTERED_DATA,'target_syllables');
                    ind=length(FILTERED_DATA.target_syllables);
                    FILTERED_DATA.target_syllables{ind+1}.Baseline_FF_mean=Baseline_FF_mean;
                    FILTERED_DATA.target_syllables{ind+1}.WN_start_daybins=WN_start_daybins;
                    FILTERED_DATA.target_syllables{ind+1}.WN_end_daybins=WNend_daybins;
                    FILTERED_DATA.target_syllables{ind+1}.syl_name=syl;
                    FILTERED_DATA.target_syllables{ind+1}.bird_num=i;
                    FILTERED_DATA.target_syllables{ind+1}.expt_num=ii;
                    FILTERED_DATA.target_syllables{ind+1}.FF_minusbase_overdays=FF_minusbase_overdays;
                    FILTERED_DATA.target_syllables{ind+1}.WN_on_day=WN_on_day;
                    FILTERED_DATA.target_syllables{ind+1}.WN_off_day=WN_off_day;
                    FILTERED_DATA.target_syllables{ind+1}.Learning_rel_targ=Learning_rel_targ;
                    FILTERED_DATA.target_syllables{ind+1}.CONSOLIDATION=CONSOLIDATION;
                    
                    FILTERED_DATA.target_syllables{ind+1}.LearnMetric_RelTarg=LearnMetric_RelTarg;
                    FILTERED_DATA.target_syllables{ind+1}.LearnMetric=LearnMetric;
                    
                    
                else
                    FILTERED_DATA.target_syllables{1}.Baseline_FF_mean=Baseline_FF_mean;
                    FILTERED_DATA.target_syllables{1}.WN_start_daybins=WN_start_daybins;
                    FILTERED_DATA.target_syllables{1}.WN_end_daybins=WNend_daybins;
                    FILTERED_DATA.target_syllables{1}.syl_name=syl;
                    FILTERED_DATA.target_syllables{1}.bird_num=i;
                    FILTERED_DATA.target_syllables{1}.expt_num=ii;
                    FILTERED_DATA.target_syllables{1}.FF_minusbase_overdays=FF_minusbase_overdays;
                    FILTERED_DATA.target_syllables{1}.WN_on_day=WN_on_day;
                    FILTERED_DATA.target_syllables{1}.WN_off_day=WN_off_day;
                    FILTERED_DATA.target_syllables{1}.Learning_rel_targ=Learning_rel_targ;
                    FILTERED_DATA.target_syllables{1}.CONSOLIDATION=CONSOLIDATION;

                    FILTERED_DATA.target_syllables{1}.LearnMetric_RelTarg=LearnMetric_RelTarg;
                    FILTERED_DATA.target_syllables{1}.LearnMetric=LearnMetric;

                end
                
                if skip_target==1;
                    % then do not put into field with all other
                    % non-targets.
                    
                    continue
                end
            end
            
            % Collect non-targets (and potentially target)
                        FILTERED_DATA.groups{group_num}.groupID=sylID_filter{2}{group_num};

            if isfield(FILTERED_DATA.groups{group_num},'syllable');
                ind=length(FILTERED_DATA.groups{group_num}.syllable);
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.Baseline_FF_mean=Baseline_FF_mean;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.WN_start_daybins=WN_start_daybins;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.WN_end_daybins=WNend_daybins;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.syl_name=syl;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.bird_num=i;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.expt_num=ii;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.FF_minusbase_overdays=FF_minusbase_overdays;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.WN_on_day=WN_on_day;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.WN_off_day=WN_off_day;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.Learning_rel_targ=Learning_rel_targ;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.CONSOLIDATION=CONSOLIDATION;
                
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.LearnMetric_RelTarg=LearnMetric_RelTarg;
                FILTERED_DATA.groups{group_num}.syllable{ind+1}.LearnMetric=LearnMetric;
                
            else
                FILTERED_DATA.groups{group_num}.syllable{1}.Baseline_FF_mean=Baseline_FF_mean;
                FILTERED_DATA.groups{group_num}.syllable{1}.WN_start_daybins=WN_start_daybins;
                FILTERED_DATA.groups{group_num}.syllable{1}.WN_end_daybins=WNend_daybins;
                FILTERED_DATA.groups{group_num}.syllable{1}.syl_name=syl;
                FILTERED_DATA.groups{group_num}.syllable{1}.bird_num=i;
                FILTERED_DATA.groups{group_num}.syllable{1}.expt_num=ii;
                FILTERED_DATA.groups{group_num}.syllable{1}.FF_minusbase_overdays=FF_minusbase_overdays;
                FILTERED_DATA.groups{group_num}.syllable{1}.WN_on_day=WN_on_day;
                FILTERED_DATA.groups{group_num}.syllable{1}.WN_off_day=WN_off_day;
                FILTERED_DATA.groups{group_num}.syllable{1}.Learning_rel_targ=Learning_rel_targ;
                FILTERED_DATA.groups{group_num}.syllable{1}.CONSOLIDATION=CONSOLIDATION;

                FILTERED_DATA.groups{group_num}.syllable{1}.LearnMetric_RelTarg=LearnMetric_RelTarg;
                FILTERED_DATA.groups{group_num}.syllable{1}.LearnMetric=LearnMetric;

            end            
        end
    end
end



%% EXTRACT SUMMARY VARIABLES

NumGroups=length(FILTERED_DATA.groups);
plotcols=lt_make_plot_colors(NumGroups,0,0);


%% FINAL LEARNING - PLOT DISTRIBUTIONS
% plot histogram of final learning values, sorted by groups (and combined)

BinEdges=[-1:0.1:1];
BinSize=BinEdges(2)-BinEdges(1);

BinCenters=BinEdges(1:end-1)+BinSize/2;


lt_figure; hold on;

% initiate data holders
Y_tot=cell(NumGroups,1); % collect data to plot

% run
% ===== SUBPLOT 1 - groups separated
lt_subplot(1,2,1); hold on;
% annotate
ylabel('count');
xlabel('Learning (rel target)');

N_string={};
for i=1:NumGroups;
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    Ncount=0;
    for ii=1:NumSylDataPts;
        
        % -- get datapoint - learning relative to target.
        Learning_rel_targ=FILTERED_DATA.groups{i}.syllable{ii}.Learning_rel_targ;
        
        
        % -- COLLECT DATA across samples (or skip if fails the filter)
        % get expt details (to decide if to collect this data
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
        issimilar=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ; % 1 if similar, 0 if not
        if extra_filter==1 && issimilar~=similar_only; % skip this if filter is on, and fails filter
            continue
        else
            
            Y_tot{i}=[Y_tot{i} Learning_rel_targ];
            
            % counter, for N
            Ncount=Ncount+1;
        end
        
    end
    
    
    
    % == Plot histogram
    [HistVals, ~]=hist(Y_tot{i}, BinCenters);
    
    
    % Ask is this population median is off from 0
    p=signrank(Y_tot{i});
    
    % plot sample size and p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+7*i, {['group ' num2str(FILTERED_DATA.groups{i}.groupID) ': N=' num2str(Ncount)],['p(median vs. 0)=' num2str(p)]}, 'FontWeight' ,'bold','FontSize',12, 'Color', plotcols{i})
    
    % normalize to get Pr
%     HistVals_Pr=HistVals/sum(HistVals);
    
    % plot
%     harea(i)=lt_plot_area(BinCenters,HistVals_Pr,plotcols{i});
%     N_string{1}=['N=' num2str(length(Y_tot{i}))]; % for sample size
    hbar(i)=lt_plot_bar(BinCenters-0.02+0.015*i, HistVals, {'Color', plotcols{i}});
    n=length(BinCenters);

    % Put line for mean
    Ymean=mean(Y_tot{i});
%     Ysem=lt_sem(Y_tot{i});    
    plot(Ymean,0,'^','Color','k','MarkerFaceColor',plotcols{i},'MarkerSize',10);
    
end

% ====== SUBPLOT 2 - all combined
lt_subplot(1,2,2); hold on;
title('All combined');
% annotate
ylabel('count');
xlabel('Learning (rel target)');

Yall=[];
for i=1:NumGroups;
    Yall=[Yall Y_tot{i}];
end

    % == Plot histogram
    [HistVals, ~]=hist(Yall, BinCenters);
    
    
    % Ask is this population median is off from 0
    p=signrank(Yall);
    
    % plot sample size and p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+7*i, {['N=' num2str(length(Yall))],['p(median vs. 0)=' num2str(p)]}, 'FontWeight' ,'bold','FontSize',12, 'Color', 'k')
    
    % normalize to get Pr
%     HistVals_Pr=HistVals/sum(HistVals);
    
    % plot
%     harea(i)=lt_plot_area(BinCenters,HistVals_Pr,plotcols{i});
%     N_string{1}=['N=' num2str(length(Y_tot{i}))]; % for sample size
    hbar(i)=lt_plot_bar(BinCenters-0.02+0.015*i, HistVals, {'Color', 'k'});
    n=length(BinCenters);

    % Put line for mean
    Ymean=mean(Yall);
%     Ysem=lt_sem(Y_tot{i});    
    plot(Ymean,0,'^','Color','k','MarkerFaceColor','k','MarkerSize',10);
    

% legend(hbar(i),{'0','1'});
lt_subtitle(['Histogram of generalization (consol start FF); data sorted by ' sylID_filter{1}]);




%% PLOT GENERALIZATION (FINAL DISTRIBUTIONS, USING LEARNING METRIC)

BinEdges=[-1:0.1:1];
BinSize=BinEdges(2)-BinEdges(1);

BinCenters=BinEdges(1:end-1)+BinSize/2;

lt_figure; hold on;

% initiate data holders
Y_tot=cell(NumGroups,1); % collect data to plot

% run
% ===== SUBPLOT 1 - groups separated
lt_subplot(1,2,1); hold on;
% annotate
ylabel('count');
xlabel('Generalization (rel target)');

N_string={};
for i=1:NumGroups;
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    Ncount=0;
    for ii=1:NumSylDataPts;
        
        % -- get datapoint - learning relative to target.
        Learning_rel_targ=FILTERED_DATA.groups{i}.syllable{ii}.LearnMetric_RelTarg;
        
        
        % -- COLLECT DATA across samples (or skip if fails the filter)
        % get expt details (to decide if to collect this data
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
        issimilar=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ; % 1 if similar, 0 if not
        if extra_filter==1 && issimilar~=similar_only; % skip this if filter is on, and fails filter
            continue
        else
            
            Y_tot{i}=[Y_tot{i} Learning_rel_targ];
            
            % counter, for N
            Ncount=Ncount+1;
        end
        
    end
    
    
    
    % == Plot histogram
    [HistVals, ~]=hist(Y_tot{i}, BinCenters);
    
    
    % Ask is this population median is off from 0
    p=signrank(Y_tot{i});
    
    % plot sample size and p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+2*i, {['group ' num2str(FILTERED_DATA.groups{i}.groupID) ': N=' num2str(Ncount)],['p(median vs. 0)=' num2str(p)]}, 'FontWeight' ,'bold','FontSize',12, 'Color', plotcols{i})
    
    % normalize to get Pr
%     HistVals_Pr=HistVals/sum(HistVals);
    
    % plot
%     harea(i)=lt_plot_area(BinCenters,HistVals_Pr,plotcols{i});
%     N_string{1}=['N=' num2str(length(Y_tot{i}))]; % for sample size
    hbar(i)=lt_plot_bar(BinCenters-0.02+0.015*i, HistVals, {'Color', plotcols{i}});
    n=length(BinCenters);

    % Put line for mean
    Ymean=mean(Y_tot{i});
%     Ysem=lt_sem(Y_tot{i});    
    plot(Ymean,0,'^','Color','k','MarkerFaceColor',plotcols{i},'MarkerSize',10);
    
    
    
end
     % ask if they are different from each other
    p=ranksum( Y_tot{1}, Y_tot{2});
    % plot p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+7, ['p(medians)=' num2str(p)], 'FontWeight' ,'bold','FontSize',12, 'Color', 'k')

% ====== SUBPLOT 2 - all combined
lt_subplot(1,2,2); hold on;
title('All combined');
% annotate
ylabel('count');
xlabel('Generalization');

Yall=[];
for i=1:NumGroups;
    Yall=[Yall Y_tot{i}];
end

    % == Plot histogram
    [HistVals, ~]=hist(Yall, BinCenters);
    
    
    % Ask is this population median is off from 0
    p=signrank(Yall);
    
    % plot sample size and p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+7*i, {['N=' num2str(length(Yall))],['p(median vs. 0)=' num2str(p)]}, 'FontWeight' ,'bold','FontSize',12, 'Color', 'k')
    
    % normalize to get Pr
%     HistVals_Pr=HistVals/sum(HistVals);
    
    % plot
%     harea(i)=lt_plot_area(BinCenters,HistVals_Pr,plotcols{i});
%     N_string{1}=['N=' num2str(length(Y_tot{i}))]; % for sample size
    hbar(i)=lt_plot_bar(BinCenters-0.02+0.015*i, HistVals, {'Color', 'k'});
    n=length(BinCenters);

    % Put line for mean
    Ymean=mean(Yall);
%     Ysem=lt_sem(Y_tot{i});    
    plot(Ymean,0,'^','Color','k','MarkerFaceColor','k','MarkerSize',10);
    
    
    

% legend(hbar(i),{'0','1'});
lt_subtitle(['Histogram of generalization, using learning metric (' Params.global.learning_metric '); data sorted by ' sylID_filter{1}]);


%% PLOT LEARNING (NOT GENERALIZATION) (FINAL DISTRIBUTIONS, USING LEARNING METRIC)

% BinEdges=[-1:0.1:1];
% BinSize=BinEdges(2)-BinEdges(1);
% 

if Params.global.learning_metric=='zscore';
    BinCenters=-3:0.3:3;
end
% BinCenters=BinEdges(1:end-1)+BinSize/2;

lt_figure; hold on;

% initiate data holders
Y_tot=cell(NumGroups,1); % collect data to plot

% run
% ===== SUBPLOT 1 - groups separated
lt_subplot(1,2,1); hold on;
% annotate
ylabel('count');
xlabel('Learning (not generalization)');

N_string={};
for i=1:NumGroups;
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    Ncount=0;
    for ii=1:NumSylDataPts;
        
        % -- get datapoint - learning relative to target.
        Learning=FILTERED_DATA.groups{i}.syllable{ii}.LearnMetric;
        
        
        % -- COLLECT DATA across samples (or skip if fails the filter)
        % get expt details (to decide if to collect this data
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
        issimilar=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ; % 1 if similar, 0 if not
        if extra_filter==1 && issimilar~=similar_only; % skip this if filter is on, and fails filter
            continue
        else
            
            Y_tot{i}=[Y_tot{i} Learning];
            
            % counter, for N
            Ncount=Ncount+1;
        end
        
    end
    
    
    
    % == Plot histogram
    [HistVals, ~]=hist(Y_tot{i}, BinCenters);
    
    
    % Ask is this population median is off from 0
    p=signrank(Y_tot{i});
    
    % plot sample size and p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+3*i, {['group ' num2str(FILTERED_DATA.groups{i}.groupID) ': N=' num2str(Ncount)],['p(median vs. 0)=' num2str(p)]}, 'FontWeight' ,'bold','FontSize',12, 'Color', plotcols{i})
    
    % normalize to get Pr
%     HistVals_Pr=HistVals/sum(HistVals);
    
    % plot
%     harea(i)=lt_plot_area(BinCenters,HistVals_Pr,plotcols{i});
%     N_string{1}=['N=' num2str(length(Y_tot{i}))]; % for sample size
    hbar(i)=lt_plot_bar(BinCenters-0.02+0.015*i, HistVals, {'Color', plotcols{i}});
    n=length(BinCenters);

    % Put line for mean
    Ymean=mean(Y_tot{i});
%     Ysem=lt_sem(Y_tot{i});    
    plot(Ymean,0,'^','Color','k','MarkerFaceColor',plotcols{i},'MarkerSize',10);
    
end

% === ADD TARGET SYLS
numsyls=length(FILTERED_DATA.target_syllables);
Yall=[];
for i=1:numsyls;
    Yall=[Yall FILTERED_DATA.target_syllables{i}.LearnMetric];
end

HistVals=hist(Yall, BinCenters);

lt_plot_bar(BinCenters, HistVals, {'Color', 'y'});

% Put line for mean
    Ymean=mean(Yall);
    plot(Ymean,0,'^','Color','k','MarkerFaceColor','y','MarkerSize',10);


% ====== SUBPLOT 2 - all combined
lt_subplot(1,2,2); hold on;
title('All combined');
% annotate
ylabel('count');
xlabel('Learning (not generalization)');

Yall=[];
for i=1:NumGroups;
    Yall=[Yall Y_tot{i}];
end

    % == Plot histogram
    [HistVals, ~]=hist(Yall, BinCenters);
    
    
    % Ask is this population median is off from 0
    p=signrank(Yall);
    
    % plot sample size and p
    Ylim=ylim;
    text(BinCenters(2), Ylim(1)+3*i, {['N=' num2str(length(Yall))],['p(median vs. 0)=' num2str(p)]}, 'FontWeight' ,'bold','FontSize',12, 'Color', 'k')
    
    % normalize to get Pr
%     HistVals_Pr=HistVals/sum(HistVals);
    
    % plot
%     harea(i)=lt_plot_area(BinCenters,HistVals_Pr,plotcols{i});
%     N_string{1}=['N=' num2str(length(Y_tot{i}))]; % for sample size
    lt_plot_bar(BinCenters, HistVals, {'Color', 'k'});
    n=length(BinCenters);

    % Put line for mean
    Ymean=mean(Yall);
%     Ysem=lt_sem(Y_tot{i});    
    plot(Ymean,0,'^','Color','k','MarkerFaceColor','k','MarkerSize',10);
    

% legend(hbar(i),{'0','1'});
lt_subtitle(['Histogram of learning, using learning metric (' Params.global.learning_metric '); data sorted by ' sylID_filter{1}]);





%% --- PLOT final learning (scatter) (WITH TEXT LABELS)
% 1) Plot final learning (each syl (per experiment) contributes one datapoint);

lt_figure; hold on;
title(['Syls sorted by: ' sylID_filter{1} '; datapoint = one syl per expt (text: bird expt targ nontarget)']);

% initiate data holders
Y_tot=cell(NumGroups,1); % collect data to plot
SylNames_tot=cell(NumGroups,1);
TargNames_tot=cell(NumGroups,1);
BirdNames_tot=cell(NumGroups,1);
ExptNames_tot=cell(NumGroups,1);

% run
for i=1:NumGroups;
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    for ii=1:NumSylDataPts;
        
        % get expt details
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
        
        % -- get data - learning relative to target.
        Learning_rel_targ=FILTERED_DATA.groups{i}.syllable{ii}.Learning_rel_targ;
        
        
        %         % troubleshooting
        %         if Learning_rel_targ>0.8;
        %             keyboard
        %         end
        
        % -- COLLECT DATA across samples (or skip if fails the filter)
        if extra_filter==1 && SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ...
                ~=similar_only; % skip this if fails filter
            continue
        else
            
            %  Collect learning
            Y_tot{i}=[Y_tot{i} Learning_rel_targ];
            
            % collect syl name for labeling purposes
            ind=length(SylNames_tot{i});
            SylNames_tot{i}{ind+1}=sylname;
            
            % collect target name
            targname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
            TargNames_tot{i}{ind+1}=targname;
            
            % collect bird name
            birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
            BirdNames_tot{i}{ind+1}=birdname;
            
            % collect expt name
            ExptID=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
            ExptNames_tot{i}{ind+1}=ExptID;
        end
    end


    % Plot all data
    lt_plot(i-0.2+0.4*rand(length(Y_tot{i}),1),Y_tot{i});
    
    xlim([0 NumGroups+1.5]);
    
    % Plot syl names
    if (1);
        for j=1:length(Y_tot{i});
            text(i+0.2,Y_tot{i}(j), [BirdNames_tot{i}{j} ' ' [ExptNames_tot{i}{j}(1:3) ExptNames_tot{i}{j}(end-2:end)]  ' ' SylNames_tot{i}{j} ' ' TargNames_tot{i}{j}], 'Color', 'b');
        end
    end
    
    % plot mean for this group
    Y_tot_mean=mean(Y_tot{i});
    Y_tot_sem=lt_sem(Y_tot{i});
    
    errorbar(i,Y_tot_mean,Y_tot_sem,'s','MarkerFaceColor','b','MarkerSize',12);
    
    lt_plot_zeroline;    
end

set(gca,'XTick',1:NumGroups);


% -- PUT X LABELS
switch sylID_filter{1}
    case 'similar_to_targ'
        set(gca,'XTickLabel',{'Different','Similar'});
    case 'preceding_syl'
        set(gca,'XTickLabel',sylID_filter{2});
end

ylabel('Learning (relative to target)');

% Perform t-test
if length(Y_tot)==2; % perform if there are two groups
    [~,p]=ttest2(Y_tot{1},Y_tot{2});
    % plot p-value
    Ylim=ylim;
    text(1.5,Ylim(2)-0.1,['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
    if p<0.05;
        plot(1.5,Ylim(2)-0.15,'*','MarkerSize',7,'Color','r');
    end
end
    

%% --- PLOT final learning (scatter) (USING LEARNING METRIC)
% 1) Plot final learning (each syl (per experiment) contributes one datapoint);

lt_figure; hold on;
title(['USing learning metric, Syls sorted by: ' sylID_filter{1} '; datapoint = one syl per expt (text: bird expt targ nontarget)']);

% initiate data holders
Y_tot=cell(NumGroups,1); % collect data to plot
SylNames_tot=cell(NumGroups,1);
TargNames_tot=cell(NumGroups,1);
BirdNames_tot=cell(NumGroups,1);
ExptNames_tot=cell(NumGroups,1);

% run
for i=1:NumGroups;
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    for ii=1:NumSylDataPts;
        
        % get expt details
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
        
        % -- get data - learning relative to target.
        Learning_rel_targ=FILTERED_DATA.groups{i}.syllable{ii}.LearnMetric_RelTarg;
        
        
        %         % troubleshooting
        %         if Learning_rel_targ>0.8;
        %             keyboard
        %         end
        
        % -- COLLECT DATA across samples (or skip if fails the filter)
        if extra_filter==1 && SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ...
                ~=similar_only; % skip this if fails filter
            continue
        else
            
            %  Collect learning
            Y_tot{i}=[Y_tot{i} Learning_rel_targ];
            
            % collect syl name for labeling purposes
            ind=length(SylNames_tot{i});
            SylNames_tot{i}{ind+1}=sylname;
            
            % collect target name
            targname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
            TargNames_tot{i}{ind+1}=targname;
            
            % collect bird name
            birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
            BirdNames_tot{i}{ind+1}=birdname;
            
            % collect expt name
            ExptID=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
            ExptNames_tot{i}{ind+1}=ExptID;
        end
    end


    % Plot all data
    lt_plot(i-0.2+0.4*rand(length(Y_tot{i}),1),Y_tot{i});
    
    xlim([0 NumGroups+1.5]);
    
    % Plot syl names
    if (1);
        for j=1:length(Y_tot{i});
            text(i+0.2,Y_tot{i}(j), [BirdNames_tot{i}{j} ' ' [ExptNames_tot{i}{j}(1:3) ExptNames_tot{i}{j}(end-2:end)]  ' ' SylNames_tot{i}{j} ' ' TargNames_tot{i}{j}], 'Color', 'b');
        end
    end
    
    % plot mean for this group
    Y_tot_mean=mean(Y_tot{i});
    Y_tot_sem=lt_sem(Y_tot{i});
    
    errorbar(i,Y_tot_mean,Y_tot_sem,'s','MarkerFaceColor','b','MarkerSize',12);
    
    lt_plot_zeroline;    
end

set(gca,'XTick',1:NumGroups);


% -- PUT X LABELS
switch sylID_filter{1}
    case 'similar_to_targ'
        set(gca,'XTickLabel',{'Different','Similar'});
    case 'preceding_syl'
        set(gca,'XTickLabel',sylID_filter{2});
end

ylabel('Learning (relative to target)');

% Perform t-test
if length(Y_tot)==2; % perform if there are two groups
    [~,p]=ttest2(Y_tot{1},Y_tot{2});
    % plot p-value
    Ylim=ylim;
    text(1.5,Ylim(2)-0.1,['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
    if p<0.05;
        plot(1.5,Ylim(2)-0.15,'*','MarkerSize',7,'Color','r');
    end
end
    

%% -- PLOT LEARNING ACROSS DAYS

% - extract stuff from params
FirstDayToPlot=Params.FilterPlot.FirstDayToPlot; % lock to what day? (1st WN day would be 1);
LastDayToPlot=Params.FilterPlot.LastDayToPlot; % lock end to what day? (5 would be 5th WN day, 'end' would be end of WN, [] would be openended);
pre_days_to_plot=3; % dont change this.


hfig=[];
hfig(1)=lt_figure; hold on; % for all data
lt_plot_format;

NumGroups=length(FILTERED_DATA.groups);
Y_tot=cell(NumGroups,1); % collect data to plot
SylNames_tot=cell(NumGroups,1);
TargNames_tot=cell(NumGroups,1);


plotcols=lt_make_plot_colors(NumGroups,0,0);

Ytot={};
Ytot_targ={};
for i=1:NumGroups;
    subplot(1,NumGroups,i); hold on;
    lt_plot_format;
    xlabel('Days');
    ylabel('Learning (ratio of learning by target');
    title(['Group:' num2str(sylID_filter{2}{i}) ]);
    
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    
    WNoffInds{i}=[]; % collect to plot
    
    for ii=1:NumSylDataPts;
        Ytot{i}{ii}=[];
        
        % collect array, one FF mean val per day
        Y=FILTERED_DATA.groups{i}.syllable{ii}.FF_minusbase_overdays; % raw, all days
        
        % make ind 4 the desired first day, taking 3 pre
        ind4=FILTERED_DATA.groups{i}.syllable{ii}.WN_on_day+FirstDayToPlot-1; % this is the day that we lock on
        indfirst=ind4-pre_days_to_plot; % is 3 days before 1st desired day
        if isempty(LastDayToPlot);
            % take last WN day
            indlast=FILTERED_DATA.groups{i}.syllable{ii}.WN_off_day;
        else
            indlast=FILTERED_DATA.groups{i}.syllable{ii}.WN_on_day+LastDayToPlot-1;
        end
        
        
        
        % -- what is learning relative to learning by target (learning last WN
        % day)
        % target info:
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
        
        FF_targ=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Learning_by_target.WN_end_daybins;
        
        % compare to target
        Y=Y./FF_targ;
        
        
        % -- COLLECT DATA & PLOT
        if extra_filter==1 && SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ...
                ~=similar_only; % skip this if fails filter
            continue
        else
            if indfirst<1;
                % then does not have enough pre days. will still lock days like
                % other data.
                Ytot{i}{ii}(1:(pre_days_to_plot-2-indfirst))=nan;
                Ytot{i}{ii}=[Ytot{i}{ii}'; Y(1:indlast)];
            else
                Ytot{i}{ii}=Y(indfirst:indlast);
            end
        end
        % collect syl name for labeling purposes
        ind=length(SylNames_tot{i});
        SylNames_tot{i}{ind+1}=sylname;
        
        % collect target name
        targname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
        TargNames_tot{i}{ind+1}=targname;
        
        
        % -- collect WN off indices for plotting
        WNoffInds{i}=[WNoffInds{i} FILTERED_DATA.groups{i}.syllable{ii}.WN_off_day-FILTERED_DATA.groups{i}.syllable{ii}.WN_on_day+1+pre_days_to_plot];
        
        % ------ PLOT this rendition
        plot(Ytot{i}{ii},'-','Color',plotcols{i});
    end


    
    % -- PLOT MEAN FOR THIS GROUP
    % what is syllable with maximum days
    tmp=[];
    for ii=1:length(Ytot{i});
        tmp=[tmp, length(Ytot{i}{ii})];
    end
    max_days=max(tmp);
    
    % make matrix with max days (dim1 = days, dim2 = renditions)
    Y_matrix=nan(max_days,length(Ytot{i}));
    for ii=1:length(Ytot{i});
        % fill in the matrix
        Y_matrix(1:length(Ytot{i}{ii}),ii)=Ytot{i}{ii};
    end
    
    % -- get mean using that matrix;
    Y_matrix_mean.nontargets{i}=nanmean(Y_matrix,2);
    Y_matrix_std.nontargets{i}=nanstd(Y_matrix,0,2);
    Y_matrix_N.nontargets{i}=sum(~isnan(Y_matrix),2);
    Y_matrix_sem.nontargets{i}=Y_matrix_std.nontargets{i}./real(sqrt(Y_matrix_N.nontargets{i}-1));
    
    % -- PLOT ALL
    shadedErrorBar(1:length(Y_matrix_mean.nontargets{i}),Y_matrix_mean.nontargets{i},Y_matrix_sem.nontargets{i},{'Color','k','LineWidth',2},1);
    
    % - Put lines to mark WN on and off
    line([pre_days_to_plot+0.5 pre_days_to_plot+0.5],ylim,'Color','b','LineWidth',2);
    % plot distribution of on and offs
    plot(WNoffInds{i},-0.7,'^b','MarkerFaceColor','b','MarkerSize',8);
    lt_plot_zeroline;
    
    
    
    % -- PLOT TRAJECTORIES FOR ALL THE TARGETS ------------------------
    NumSylDataPts=length(FILTERED_DATA.target_syllables);
    for ii=1:NumSylDataPts;
        Ytot_targ{ii}=[];
        
        % collect array, one FF mean val per day
        Y=FILTERED_DATA.target_syllables{ii}.FF_minusbase_overdays; % raw, all days
        
        % make ind 4 the desired first day, taking 3 pre
        ind4=FILTERED_DATA.target_syllables{ii}.WN_on_day+FirstDayToPlot-1; % this is the day that we lock on
        indfirst=ind4-pre_days_to_plot; % is 3 days before 1st desired day
        if isempty(LastDayToPlot);
            % take last WN day
            indlast=FILTERED_DATA.target_syllables{ii}.WN_off_day;
        else
            indlast=FILTERED_DATA.target_syllables{ii}.WN_on_day+LastDayToPlot-1;
        end
        
        % -- what is learning relative to learning by target (learning last WN
        % day)
        % target info:
        birdnum=FILTERED_DATA.target_syllables{ii}.bird_num;
        exptnum=FILTERED_DATA.target_syllables{ii}.expt_num;
        sylname=FILTERED_DATA.target_syllables{ii}.syl_name;
        
        FF_targ=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Learning_by_target.WN_end_daybins;
        
        % compare to target
        Y=Y./FF_targ;
        
        % -- COLLECT DATA
            if indfirst<1;
                % then does not have enough pre days. will still lock days like
                % other data.
                Ytot_targ{ii}(1:(pre_days_to_plot-2-indfirst))=nan;
                Ytot_targ{ii}=[Ytot_targ{ii}'; Y(1:indlast)];
            else
                Ytot_targ{ii}=Y(indfirst:indlast);
            end
        
        % ------ PLOT this rendition
        plot(Ytot_targ{ii},'-','Color','b');
    end

    % - Plot mean
    % what is syllable with maximum days
    tmp=[];
    for ii=1:length(Ytot_targ);
        tmp=[tmp, length(Ytot_targ{ii})];
    end
    max_days=max(tmp);
    
    % make matrix with max days (dim1 = days, dim2 = renditions)
    Y_matrix=nan(max_days,length(Ytot_targ));
    for ii=1:length(Ytot_targ);
        % fill in the matrix
        Y_matrix(1:length(Ytot_targ{ii}),ii)=Ytot_targ{ii};
    end
    
    % -- get mean using that matrix;
    Y_matrix_mean.targets=nanmean(Y_matrix,2);
    Y_matrix_std.targets=nanstd(Y_matrix,0,2);
    Y_matrix_N.targets=sum(~isnan(Y_matrix),2);
    Y_matrix_sem.targets=Y_matrix_std.targets./real(sqrt(Y_matrix_N.targets-1));
    
    % -- plot
    shadedErrorBar(1:length(Y_matrix_mean.targets),Y_matrix_mean.targets,Y_matrix_sem.targets,{'Color','b','LineWidth',2},1);

    % -----------------------------------------------------------------
    
end

lt_subtitle(['[CONSOL EARLY FF] Sorted by: ' sylID_filter{1} ': each syl (per expt) one line']);


% -- SUMMARY PLOTS
lt_figure; hold on;

% -- nontargets
for i=1:NumGroups;
    lt_subplot(1,NumGroups,i); hold on;
    xlabel('Days');
    ylabel('Learning (ratio of learning by target');
    
    title(['Group:' num2str(sylID_filter{2}{i}) ]);
    % -- plot
    shadedErrorBar(1:length(Y_matrix_mean.nontargets{i}),Y_matrix_mean.nontargets{i},Y_matrix_sem.nontargets{i},{'Color',plotcols{i},'LineWidth',2},1);
    
    % -- target
    shadedErrorBar(1:length(Y_matrix_mean.targets),Y_matrix_mean.targets,Y_matrix_sem.targets,{'Color','b','LineWidth',2},1);
    
    % plot distribution of on and offs
    plot(WNoffInds{i},-0.7,'^b','MarkerFaceColor','b','MarkerSize',8);
    % - Put lines to mark WN on and off
    line([pre_days_to_plot+0.5 pre_days_to_plot+0.5],ylim,'Color','b','LineWidth',2);
    lt_plot_zeroline;
end
lt_subtitle(['[CONSOL EARLY FF] Sorted by: ' sylID_filter{1} ': each syl (per expt) one line']);



%% === CONSOLIDATION - plot all data 

% figure notification
lt_figure;
lt_plot_text(0.1,0.1, ['NOTE: CONSOLIDATION IS ALL USING EARLY CONSOL FF, NOT LEARNING METRIC'])

% ==== CONTINUE
MinConsolidDur=Params.FilterPlot.MinConsolidDur;
MaxDayFromWNStart=Params.FilterPlot.MaxDayFromWNStart; 


disp(' ');
% == SCATTER OF CHANGE IN LEARNING VS. EARLY LEARNING.
Xtot={};
Ytot={};
for i=1:NumGroups;
    NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
    
    % -- COLLECT DATA
    Xtot{i}=[];
    Ytot{i}=[];
    for ii=1:NumSylDataPts;
        
        % -------- THROW OUT SOME DATA
        % Only keep datapoint from experiments with sufficiently long
        % consolidation periods
        birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
        exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
        birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
        exptID=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
        
        ConsolidDur=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.ConsolEndInd - ...
            SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.ConsolStartInd+1;
        
        if ConsolidDur<MinConsolidDur;
            disp(['Skipped group ' num2str(i) ' datapoint: ' num2str(ii) ' (bird ' birdname ', experiment ' exptID '); lacking enough consolidation days']);
            continue
        end
        
        % Only keep datapoint if consolidation started early enough after
        % WN start
        ConsolStartInd=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.ConsolStartInd;
        WNStartInd=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        
        DaysWNtoConsolid=ConsolStartInd-WNStartInd+1;
        
        if DaysWNtoConsolid>MaxDayFromWNStart;
            disp(['Skipped group ' num2str(i) ' datapoint: ' num2str(ii) ' (bird ' birdname ', experiment ' exptID '); consolidation started too late']);
            %             keyboard
            continue
        end
        
        
        % ---------- COLLECT DATA
        % early learning
        X=FILTERED_DATA.groups{i}.syllable{ii}.CONSOLIDATION.Learning_start; % early learning
        
        % change in learning
        Y=FILTERED_DATA.groups{i}.syllable{ii}.CONSOLIDATION.Learning_end-FILTERED_DATA.groups{i}.syllable{ii}.CONSOLIDATION.Learning_start; % late minus early
        
        % collect across all syls
        Xtot{i}=[Xtot{i} X];
        Ytot{i}=[Ytot{i} Y];
    end
end

% === PLOT DOTS (early vs late)
lt_figure;
hsplot=[];
Yall=[];
for i=1:NumGroups;
    hsplot(i)=lt_subplot(1,3,i); hold on;
    title(['Group ' num2str(i)]);
    
    X=[1 2];
    Y=[Xtot{i}' Xtot{i}'+Ytot{i}'];
    
    lt_plot(X,Y,{'LineStyle','-'});
    
    xlim([-0.5 3.5]);
    
    % plot mean
    errorbar([1.2 2.2], mean(Y), [lt_sem(Y(:,1)) lt_sem(Y(:,2))],'-s','MarkerSize',10);
    
    % ttest
    [p, h]=signrank(Y(:,1),Y(:,2));
    
    ylim([-0.7 0.7])
    lt_plot_zeroline
    
    % plot p-value
    lt_plot_pvalue(p)
    
    % collect for all data
    Yall=[Yall; Y];
end

% --- FOR ALL DATA
lt_subplot(1,3,3); hold on;
title('All combined');

lt_plot(X,Yall,{'LineStyle','-'});
xlim([-0.5 3.5]);
ylim([-0.7 0.7])
lt_plot_zeroline

% plot mean
errorbar([1.2 2.2], mean(Yall), [lt_sem(Yall(:,1)) lt_sem(Yall(:,2))],'-s','MarkerSize',10);

% test signi
[p, h]=signrank(Yall(:,1),Yall(:,2));

lt_plot_pvalue(p)

% -------
lt_subtitle('Learning (early vs late)');



% == PLOT DOTS (absolute value of learning) - ASKS whether consolidation
% tends to drive syllables closer to baseline (but NOT necessarily whether
% drives in direction - e.g. can drive back to reversion and overshoot, and
% would be counted as no change in absolulte value)
lt_figure;
hsplot=[];
Yall=[];
for i=1:NumGroups;
    lt_subplot(1,3,i); hold on;
    title(['Group ' num2str(i)]);
    
    X=[1 2];
    Y=[abs(Xtot{i}') abs(Xtot{i}'+Ytot{i}')];
    
    lt_plot(X,Y,{'LineStyle','-'});
    
    xlim([-0.5 3.5]);
ylim([0 0.8])
lt_plot_zeroline    
    
    % plot mean
    errorbar([1.2 2.2], mean(Y), [lt_sem(Y(:,1)) lt_sem(Y(:,2))],'-s','MarkerSize',10);

    % ttest
    [p, h]=signrank(Y(:,1),Y(:,2));   
    

% plot p-value
lt_plot_pvalue(p)

% collect for all;
Yall=[Yall; Y];

end

% -- PLOT ALL
lt_subplot(1,3,3); hold on;
title('All combined');

lt_plot(X,Yall,{'LineStyle','-'});
    xlim([-0.5 3.5]);
ylim([0 0.8])
lt_plot_zeroline    

% plot mean
errorbar([1.2 2.2], mean(Yall), [lt_sem(Yall(:,1)) lt_sem(Yall(:,2))],'-s','MarkerSize',10);

% test signi
[p, h]=signrank(Yall(:,1),Yall(:,2));

lt_plot_pvalue(p)

% -------

lt_subtitle('Learning (absolute value), early vs late');



% == PLOT DOTS (flipping sign of data with initial negative learning (flips both early and late, if early is negative)
% - analogous to looking at whether magnitude
% of early learning predicts magnitude of reversion to baseline) - an issue
% is that overshoots will be counted as reversion, even though is a little
% different.) - ASKS whether early learning predicts magnitude of change in
% direction of baseline.

lt_figure;
Yall=[];
for i=1:NumGroups;
    lt_subplot(1,3,i); hold on;
    title(['Group ' num2str(i)]);
    ylabel('learning (sign flipped if early learning is neg)');
    
    inds_to_flip=(Xtot{i}<0);
    
    % replace early learning
    X_replace=Xtot{i}.*-inds_to_flip;
    X=Xtot{i};
    X(inds_to_flip)=X_replace(inds_to_flip);
    
    % flip late learning too
    Y_replace=Ytot{i}.*-inds_to_flip; % these are negative early learning, flip late learning of same syls
    Y=Ytot{i};
    Y(inds_to_flip)=Y_replace(inds_to_flip);
    
    Y=Y+X; % add, because Y is difference
    
    Y=[X' Y']; % [early late];
    X=[1 2];
    
    lt_plot(X,Y,{'LineStyle','-'});
    
    xlim([-0.5 3.5]);
ylim([-0.5 0.8])
lt_plot_zeroline    
    
    % plot mean
    errorbar([1.2 2.2], mean(Y), [lt_sem(Y(:,1)) lt_sem(Y(:,2))],'-s','MarkerSize',10);

    % ttest
    [p, h]=signrank(Y(:,1),Y(:,2));   
    
% plot p-value
lt_plot_pvalue(p)

% collect for all;
Yall=[Yall; Y];

end

% -- PLOT ALL
lt_subplot(1,3,3); hold on;
title('All combined');

lt_plot(X,Yall,{'LineStyle','-'});
    xlim([-0.5 3.5]);
ylim([0 0.8])
lt_plot_zeroline    

% plot mean
errorbar([1.2 2.2], mean(Yall), [lt_sem(Yall(:,1)) lt_sem(Yall(:,2))],'-s','MarkerSize',10);

% test signi
[p, h]=signrank(Yall(:,1),Yall(:,2));

lt_plot_pvalue(p)

% -------

lt_subtitle('Magnitude of change in direction of unlearning, early vs late');




% == PLOT SCATTER (compare change vs. start)
lt_figure;
hsplot=[];
Xall=[];
Yall=[];
for i=1:NumGroups;
    lt_subplot(1,3,i); hold on;
    title(['Group ' num2str(i)]);
    ylabel('Change (difference)');
    xlabel('Learning, start of consolid');
    
    X=Xtot{i}'; % early
    Y=Ytot{i}'; % late minus early
    
    lt_plot(X,Y,{'Color',plotcols{i}});
    
    % regression
    [~,~,~,~,~,regr_stats]=lt_regress(Y ,X,1,0);
    

    line([0 0], ylim);
    line(xlim, [0 0]);

    % collect for all
    Xall=[Xall; X];
    Yall=[Yall; Y];
end

% -- PLOT FOR ALL
    lt_subplot(1,3,3); hold on;
    title(['All combined']);
    ylabel('Change (difference)');
    xlabel('Learning, start of consolid');
    
    lt_plot(Xall,Yall);
    
    % regression
    [~,~,~,~,~,regr_stats]=lt_regress(Yall, Xall,1,0);
    
    line([0 0], ylim);
    line(xlim, [0 0]);

% ---

lt_subtitle(['Change over consolidation (Group: ' num2str(FILTERED_DATA.groups{i}.groupID) ')']);






% == PLOT scatter of pre vs post
lt_figure;
hsplot=[];
Xall=[];
Yall=[];
for i=1:NumGroups;
    hsplot(i)=lt_subplot(1,3,i); hold on;
    title(['Group ' num2str(i)]);
    xlabel('Early learning');
    ylabel('Late learning');
    
    X=Xtot{i}'; % early
    Y=Xtot{i}'+Ytot{i}'; % late
    
    lt_plot(X,Y,{'Color',plotcols{i}});
    

    xlim([-1 1]);
    ylim([-1 1]);
    line([0 0], ylim);
    line(xlim, [0 0]);
    line(xlim,ylim);
    
    % collect all
    Xall=[Xall; X];
    Yall=[Yall; Y];

end

% -- PLOT ALL
    lt_subplot(1,3,3); hold on;
    title('All combined');
    xlabel('Early learning');
    ylabel('Late learning');
    
    lt_plot(Xall,Yall);
    
    xlim([-1 1]);
    ylim([-1 1]);
    line([0 0], ylim);
    line(xlim, [0 0]);
    line(xlim,ylim);
    
% ---

lt_subtitle('Early vs. late learning');



% TO PLOT
% 4) Plot with same day windows

%% CONSOLIDATION - plot change, dividing groups into positive and negative generalizers
% == PLOT SCATTER (compare change vs. start)
lt_figure;
hsplot=[];
Xall=[];
Yall=[];
for i=1:NumGroups;
    lt_subplot(1,3,i); hold on;
    ylabel('Change over consolid (hz)');
    title(['group num ' num2str(i)]);
    
%     title(['Group ' num2str(i)]);
%     ylabel('Change (difference)');
%     xlabel('Learning, start of consolid');
    
    X=Xtot{i}'; % early
    Y=Ytot{i}'; % late minus early
    
    % --- Negative generalizers
    inds=X<0;
    Y_neg=Y(inds);

    % --- Positive generalizers
    inds=X>0;
    Y_pos=Y(inds);
    
    % -------- PLOT SCATTER
    lt_plot(1, Y_neg) % negative
    lt_plot(2, Y_pos) % positive
    xlim([0 3]);
   
    % ---------- PLOT MEANS
    errorbar(1.1, mean(Y_neg), lt_sem(Y_neg), 'sb', 'MarkerSize', 9) % negative
    errorbar(2.1, mean(Y_pos), lt_sem(Y_pos), 'sb', 'MarkerSize', 9) % positive
    
    % ------ significance test
    p=ranksum(Y_neg, Y_pos);
    lt_plot_pvalue(p);
    
    lt_plot_zeroline;
    set(gca, 'XTick', 1:2);
    set(gca, 'XTickLabel', {'start neg' ,'start pos'})
    
    % collect for all
    Xall=[Xall; X];
    Yall=[Yall; Y];
end

% -- PLOT FOR ALL
    lt_subplot(1,3,3); hold on;

    ylabel('Change over consolid (hz)');
    title('all');

    % --- Negative generalizers
    inds=Xall<0;
    Y_neg=Yall(inds);

    % --- Positive generalizers
    inds=Xall>0;
    Y_pos=Yall(inds);
    
    % -------- PLOT SCATTER
    lt_plot(1, Y_neg) % negative
    lt_plot(2, Y_pos) % positive
    xlim([0 3]);
   
    % ---------- PLOT MEANS
    errorbar(1.1, mean(Y_neg), lt_sem(Y_neg), 'sb', 'MarkerSize', 9) % negative
    errorbar(2.1, mean(Y_pos), lt_sem(Y_pos), 'sb', 'MarkerSize', 9) % positive
    
    % ------ significance test
    p=ranksum(Y_neg, Y_pos);
    lt_plot_pvalue(p);
    
    lt_plot_zeroline;

% ---
lt_subtitle(['Change over consolidation, serparating into neg and pos generalizers']);


%% CONSOLIDATION - plot change, dividing groups greater and less than median
% == PLOT SCATTER (compare change vs. start)
lt_figure;
hsplot=[];
Xall=[];
Yall=[];
for i=1:NumGroups;
    lt_subplot(1,3,i); hold on;
    ylabel('Change over consolid (hz)');
    title(['group num ' num2str(i)]);
    
%     title(['Group ' num2str(i)]);
%     ylabel('Change (difference)');
%     xlabel('Learning, start of consolid');
    
    X=Xtot{i}'; % early
    Y=Ytot{i}'; % late minus early
    
    % divider is median generalization score
    median_generalization=median(X);
    
    % --- Negative generalizers
    inds=X<median_generalization;
    Y_neg=Y(inds);

    % --- Positive generalizers
    inds=X>median_generalization;
    Y_pos=Y(inds);
    
    % -------- PLOT SCATTER
    lt_plot(1, Y_neg) % negative
    lt_plot(2, Y_pos) % positive
    xlim([0 3]);
   
    % ---------- PLOT MEANS
    errorbar(1.1, mean(Y_neg), lt_sem(Y_neg), 'sb', 'MarkerSize', 9) % negative
    errorbar(2.1, mean(Y_pos), lt_sem(Y_pos), 'sb', 'MarkerSize', 9) % positive
    
    % ------ significance test
    p=ranksum(Y_neg, Y_pos);
    lt_plot_pvalue(p);
    
    lt_plot_zeroline;
    set(gca, 'XTick', 1:2);
    set(gca, 'XTickLabel', {'start neg' ,'start pos'})
    
    % collect for all
    Xall=[Xall; X];
    Yall=[Yall; Y];
end

% -- PLOT FOR ALL
    lt_subplot(1,3,3); hold on;

    ylabel('Change over consolid (hz)');
    title('all');

        % divider is median generalization score
    median_generalization=median(X);
    
    % --- Negative generalizers
    inds=Xall<median_generalization;
    Y_neg=Yall(inds);

    % --- Positive generalizers
    inds=Xall>median_generalization;
    Y_pos=Yall(inds);
    
    % -------- PLOT SCATTER
    lt_plot(1, Y_neg) % negative
    lt_plot(2, Y_pos) % positive
    xlim([0 3]);
   
    % ---------- PLOT MEANS
    errorbar(1.1, mean(Y_neg), lt_sem(Y_neg), 'sb', 'MarkerSize', 9) % negative
    errorbar(2.1, mean(Y_pos), lt_sem(Y_pos), 'sb', 'MarkerSize', 9) % positive
    
    % ------ significance test
    p=ranksum(Y_neg, Y_pos);
    lt_plot_pvalue(p);
    
    lt_plot_zeroline;

% ---
lt_subtitle(['Change over consolidation, serparating into generalization < or > than median (within group)']);

%% ================================== CHANGE IN TARGET OVER CONSOLIDATION

numtargs=length(FILTERED_DATA.target_syllables);

% Collect all starts and ends
Yend=[];
Ystart=[];
    Ysem_start=[];
    Ysem_end=[];

for i=1:numtargs;
    
    Ystart=[Ystart FILTERED_DATA.target_syllables{i}.CONSOLIDATION.FFminusbase_start];
    Yend=[Yend FILTERED_DATA.target_syllables{i}.CONSOLIDATION.FFminusbase_end];
    
    Ysem_start=[Ysem_start FILTERED_DATA.target_syllables{i}.CONSOLIDATION.FFsem_start];
    Ysem_end=[Ysem_end FILTERED_DATA.target_syllables{i}.CONSOLIDATION.FFsem_end];
    
end

% ==== PLOT
lt_figure; hold on;
title('Change over consolidation by target syls');
ylabel('FF (hz minus base)');

for i=1:length(Ystart);
    errorbar([1 2], [Ystart(i) Yend(i)], [Ysem_start(i) Ysem_end(i)], 'o-', 'MarkerSize',8);
end
xlim([0 3]);

% plot mean
YYmean=mean([Ystart' Yend']);
YYsem=[lt_sem(Ystart) lt_sem(Yend)];
    errorbar([1.2 2.2], YYmean, YYsem, 's-', 'Color' ,'k', 'MarkerFaceColor','k', 'MarkerSize',8);

% p-value
p=signrank(Ystart, Yend);
lt_plot_pvalue(p)
lt_plot_zeroline;

% ==== PLOT (absolute value)
lt_figure; hold on;
title('Change over consolidation by target syls (abs value)');
ylabel('FF (hz minus base)');

for i=1:length(Ystart);
    errorbar([1 2], [abs(Ystart(i)) abs(Yend(i))], [Ysem_start(i) Ysem_end(i)], 'o-', 'MarkerSize',8);
end
xlim([0 3]);

% plot mean
YYmean=mean(abs([Ystart' Yend']));
YYsem=[lt_sem(abs(Ystart)) lt_sem(abs(Yend))];
    errorbar([1.2 2.2], YYmean, YYsem, 's-', 'Color' ,'k',  'MarkerFaceColor','k','MarkerSize',8);

% p-value
p=signrank(abs(Ystart), abs(Yend));
lt_plot_pvalue(p)


%% BELOW - REPLACED WITH CONSOLIDATION ANALYSIS...
% --- PLOT, COMPARE EARLY TO LATE LEARNING ---------------------------
% very crude, just take 1st 5 days as early, and rest as late
% figure; hold on;
% lt_plot_format;
% title(['Learning relative to target syl. Data sorted by: ' sylID_filter{1}]);
% 
% % -- for targets, 
%     % number of total days
%     TotalWNDays=length(Y_matrix_mean.targets)-pre_days_to_plot;
%     
%     % Early learning
%     earlyinds=pre_days_to_plot+1:pre_days_to_plot+5;
%     Y_TargEarly=mean(Y_matrix_mean.targets(earlyinds));
%         
%     % late learning
%     lateinds=pre_days_to_plot+6:TotalWNDays;
%     Y_TargLate=mean(Y_matrix_mean.targets(lateinds));
%     
% 
% % -- nontargets
% hplot=[];
% for i=1:NumGroups;
%     % number of total days
%     TotalWNDays=length(Y_matrix_mean.targets)-pre_days_to_plot;
%     
%     % Early learning
%     earlyinds=pre_days_to_plot+1:pre_days_to_plot+5;
%     Yratio.early{i}=mean(Y_matrix_mean.nontargets{i}(earlyinds))/Y_TargEarly;
%         
%     % late learning
%     lateinds=pre_days_to_plot+6:TotalWNDays;
%     Yratio.late{i}=mean(Y_matrix_mean.nontargets{i}(lateinds))/Y_TargLate;
%     
%     
%     % - Plot
%     X=[1 2];
%     hplot(i)=plot(X, [Yratio.early{i} Yratio.late{i}],'-o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',10);
%     xlim([0.5 2.5]);
%     
%     set(gca,'XTick', [1 2]);
%     set(gca,'XTickLabel',{'Early','Late'})
%     ylabel('Learning (deltaHz relative to target syl)');
% end
% % legend(hplot,['Group: ' numstr(sylID_filter{2})]);
% lt_plot_zeroline



