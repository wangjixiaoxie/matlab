function [StimEpochs_aligned]=lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs1(PARAMS, DATSTRUCT, StimEpochs, RunBin, NumEdgeTrials, TimeFieldsOfInterest, statfield)


%% PLOTS ALL STIM EPOCHS IN ONE PLOT (ACROSS DAYS), THEN PLOTS ALL INDIVIDUALLY
NumDays=length(DATSTRUCT.data);
NumTrials_plot=[]; % keep empty to get all.

TrialTypes={'StimCatch','StimNotCatch','All_preceding'};
% NumTrials=30; % e.g. of stim catch get 20. (make empty if want all)

    for i=1:length(TimeFieldsOfInterest);
        hfig=lt_figure; hold on; title(['All data (no Baseline), divided into stim epochs; timefield ' num2str(i)]);
%         h(2)=figure; hold on; title(['timefield ' num2str(i) 'raw values and means']);
        
        for ii=1:NumDays;
            if ~isempty(StimEpochs{ii}); % this day has stims
                
                NumEpochs=length(StimEpochs{ii}.epoch);
                for k=1:NumEpochs;
                    
                    for kk=1:length(TrialTypes)
                        trialtype=TrialTypes{kk};
                        
                        try % because some days without preceding non-stim
                            Inds=StimEpochs{ii}.epoch{k}.(trialtype);
                            
                            % take edge trials if wanted
                            if ~isempty(NumTrials_plot);
                                if length(Inds)>=NumTrials_plot; % if has enough trials, else keeps Inds as is.
                                    switch trialtype
                                        case 'All_preceding'
                                            Inds=Inds(end-NumTrials_plot+1:end);
                                        otherwise
                                            Inds=Inds(1:NumTrials_plot);
                                    end
                                end
                            end
                            
                            switch trialtype % becuase names of fields diffeernt.
                                case 'All_preceding';
                                    trialtype='All';
                            end
                            
                            
                            Yvals=DATSTRUCT.data{ii}.timewindow{i}.(trialtype).MINUSBaseHrBins.([statfield '_fudgeDST'])(Inds);
                            Tvals=DATSTRUCT.data{ii}.timewindow{i}.(trialtype).MINUSBaseHrBins.timevals_fudgeDST(Inds);
                        catch err
                        end
                        
                        % 1) PLOT RUNNING AVG
%                         figure(h(1));
                        
                        Yvals_sm=lt_running_stats(Yvals, RunBin);
                        Tvals_sm=lt_running_stats(Tvals, RunBin);
                        
                        % plot
                        switch trialtype
                            case 'StimCatch'
                                shadedErrorBar(ii+Tvals_sm.Median./24,Yvals_sm.Mean,Yvals_sm.SEM,{'-r','LineWidth',2},1);
                                %                             plot(ii+Tvals_sm.Median./24,Yvals_sm.Mean,'-r','LineWidth',2);
                            case 'StimNotCatch'
                                shadedErrorBar(ii+Tvals_sm.Median./24,Yvals_sm.Mean,Yvals_sm.SEM,{'-g','LineWidth',2},1);
                                %                             plot(ii+Tvals_sm.Median./24,Yvals_sm.Mean,'-g');
                            case 'All'
                                shadedErrorBar(ii+Tvals_sm.Median./24,Yvals_sm.Mean,Yvals_sm.SEM,{'-k','LineWidth',2},1);
                                %                             plot(ii+Tvals_sm.Median./24,Yvals_sm.Mean,'-k');
                        end
                        
                        
                        % 2) get means and plot
%                         figure(h(2));
                        
                        Yvals_mean=mean(Yvals);
                        Yvals_sem=lt_sem(Yvals);
                        
                        Tvals_median=median(Tvals);
                        
                        % plot raw vals
                        switch trialtype
                            case 'StimCatch'
                                plot(ii+Tvals./24,Yvals,'.r');
                            case 'StimNotCatch'
                                plot(ii+Tvals./24,Yvals,'.g')
                            case 'All'
                                plot(ii+Tvals./24,Yvals,'.k')
                        end
                        
                        
                        % plot means
                        switch trialtype
                            case 'StimCatch'
                                errorbar(ii+Tvals_median./24,Yvals_mean,Yvals_sem,'ok','MarkerFaceColor','r','MarkerSize',10);
                            case 'StimNotCatch'
                                errorbar(ii+Tvals_median./24,Yvals_mean,Yvals_sem,'ok','MarkerFaceColor','g','MarkerSize',10);
                            case 'All'
                                errorbar(ii+Tvals_median./24,Yvals_mean,Yvals_sem,'ok','MarkerFaceColor','k','MarkerSize',10);
                        end
                        
                    end
                end
            end
        end
        line(xlim,[0 0],'Color','k','LineStyle','--'); % baseline line
        
    end
    
    
    %% PLOT ALL STIM EPOCHS ALONE
    % === COLLECT ALL STIM EPOCHS, aligned to onset of stim.
StimEpochs_aligned={};
TrialTypes2={'StimCatch','StimNotCatch','All_preceding','All_post'};

for i=1:length(TimeFieldsOfInterest);
    c=1;
    for ii=1:NumDays;
        % remove baseline days
        if ~any(ii==PARAMS.global.BaselineDays)
            if ~isempty(StimEpochs{ii}); % day has stim epochs
                
                NumEpochs=length(StimEpochs{ii}.epoch); % how many epochs today
                
                for k=1:NumEpochs
                    
                    if isfield(StimEpochs{ii}.epoch{k},'All_preceding'); % make sure this epoch also has pre-stim data
                        
                        % COLLECT INFORMATION ABOUT THIS EPOCH
                        % 0) First get t0 (i.e. the time of the first stim
                        % catch or notcatch).
                        
                        Inds=StimEpochs{ii}.epoch{k}.StimCatch;
                        Tvals=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.MINUSBaseHrBins.timevals_fudgeDST(Inds);
                        
                        tmp1=min(Tvals);
                        
                        
                        Inds=StimEpochs{ii}.epoch{k}.StimNotCatch;
                        Tvals=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.timevals_fudgeDST(Inds);
                        
                        tmp2=min(Tvals);
                        
                        time_zero=min([tmp1 tmp2]);
                        
                        % Put info about epoch into structure
                        StimEpochs_aligned.timewindow{i}.epoch{c}.epoch_info.time_zero=time_zero;
                        StimEpochs_aligned.timewindow{i}.epoch{c}.epoch_info.day=ii;
                        StimEpochs_aligned.timewindow{i}.epoch{c}.epoch_info.withinday_epoch=k;
                        
                        
                        
                        
                        % COLLECT DATA FOR THIS EPOCH
                        for iii=1:length(TrialTypes2);
                            trialtype=TrialTypes2{iii};
                            
                            % get inds
                            if isfield(StimEpochs{ii}.epoch{k},trialtype)
                                Inds=StimEpochs{ii}.epoch{k}.(trialtype);
                                
                                % how many datapoints to take?
                                if ~isempty(NumTrials_plot);
                                    % if there are not enough data, give user warning
                                    if length(Inds)<NumTrials_plot;
                                        disp(['warning, not enough trials for ' trialtype ' for day: ' num2str(ii), ' epoch: ' num2str(k) ' have ' num2str(length(Inds)) ' trials.']);
                                    else % have enough trials
                                        % get data from end or beginning of epoch
                                        switch trialtype
                                            case 'All_preceding'
                                                Inds=Inds(end-NumTrials_plot+1:end);
                                            otherwise
                                                Inds=Inds(1:NumTrials_plot);
                                        end
                                    end
                                end
                                
                                    % might have to change name of trialtype.
                                    switch trialtype
                                        case {'All_preceding', 'All_post'}
                                            trialtype2='All';
                                        otherwise
                                            trialtype2=trialtype;
                                    end
                                
                                % raw data
                                Yvals=DATSTRUCT.data{ii}.timewindow{i}.(trialtype2).MINUSBaseHrBins.([statfield '_fudgeDST'])(Inds);
                                Tvals=DATSTRUCT.data{ii}.timewindow{i}.(trialtype2).MINUSBaseHrBins.timevals_fudgeDST(Inds);
                                
                                % convert times to times relative to time zero
                                Tvals_minus_t0=Tvals-time_zero;
                                
                                % Put into structure
                                StimEpochs_aligned.timewindow{i}.epoch{c}.data.(trialtype).([statfield '_fudgeDST'])=Yvals;
                                StimEpochs_aligned.timewindow{i}.epoch{c}.data.(trialtype).timevals_fudgeDST_minust0=Tvals_minus_t0;
                                
                            end
                            
                        end
                        % ITERATE ONE TIME
                        c=c+1;
                    end
                end
            end
        end
    end
end


% PLOT
for k=1:length(TimeFieldsOfInterest)
    
    NumStimEpochs=length(StimEpochs_aligned.timewindow{k}.epoch);
    
    % PLOT ALL STIM EPOCHS ALONE
    lt_figure; hold on;
    
    
    for i=1:NumStimEpochs;
        hsplot(i)=subplot(ceil(sqrt(NumStimEpochs)),ceil(sqrt(NumStimEpochs)),i); hold on;
        title(['Day: ' num2str(StimEpochs_aligned.timewindow{k}.epoch{i}.epoch_info.day) ', Time: ' num2str(StimEpochs_aligned.timewindow{k}.epoch{i}.epoch_info.time_zero) 'hrs']);
        
        for ii=1:length(TrialTypes2);
            trialtype=TrialTypes2{ii};
            
            if isfield(StimEpochs_aligned.timewindow{k}.epoch{i}.data,trialtype);
                
                Tvals=StimEpochs_aligned.timewindow{k}.epoch{i}.data.(trialtype).timevals_fudgeDST_minust0;
                Yvals=StimEpochs_aligned.timewindow{k}.epoch{i}.data.(trialtype).([statfield '_fudgeDST']);
                
                % perform local regression to smooth
                
                Yvals_sm=smooth(Yvals, 15, 'rlowess');
                
                
                
                % PLOT
                switch trialtype
                    case {'All_preceding', 'All_post'}
                        plotcol='k';
                    case 'StimCatch'
                        plotcol='r';
                    case 'StimNotCatch'
                        plotcol='g';
                end
                
                % x is time
                plot(Tvals,Yvals,'.','Color',plotcol)
                plot(Tvals,Yvals_sm,'-','Color',plotcol, 'LineWidth', 2)
                
                
                %             % Plot - x is rends (pre) and time (post)
                %             switch trialtype
                %                 case 'All_preceding'
                %                     %                 Xrends=-length(Tvals):-1;
                %                     Xrends=-linspace(0,max(abs(Tvals)),length(Tvals)); % take uniform numbers from negative (time of start) to 0.
                %
                %                 otherwise
                %                     Xrends=Tvals; % use time vales
                %             end
                %
                %             plot(Xrends,Yvals,'.','Color',plotcol)
                %             plot(Xrends,Yvals_sm,'-','Color',plotcol)
                
                
                % Plot - x is rends (all)
                if (0)
                    switch trialtype
                        case 'All_preceding'
                            Xrends=-length(Tvals):-1;
                        case 'All_post'
                            % how many renditions has the stim trials taken up?
                            tmp=max([length(StimEpochs_aligned.timewindow{k}.epoch{i}.data.StimCatch.timevals_fudgeDST_minust0), ...
                                length(StimEpochs_aligned.timewindow{k}.epoch{i}.data.StimNotCatch.timevals_fudgeDST_minust0)]);
                            Xrends=tmp+1:tmp+length(Tvals);
                        otherwise
                            Xrends=1:length(Tvals); % use rends
                            
                    end
                    
                    plot(Xrends,Yvals,'.','Color',plotcol)
                    plot(Xrends,Yvals_sm,'-','Color',plotcol)
                end
            end
        end
    end
    
    linkaxes(hsplot,'xy')
    ylim([-500 500]);
    lt_subtitle(['timewindow ' num2str(k)]);
end



%% Are stim catch trials accurate representations of non-stim trials?

NumTrials2=NumEdgeTrials; % average over N trials at edges

% for each stim epoch, get four numbers (2 per transition).
for i=1:length(TimeFieldsOfInterest);
    
    NumEpochs=length(StimEpochs_aligned.timewindow{i}.epoch); % how many epochs today
    
    for k=1:NumEpochs
        
        
        
        % extract raw values into matrix
        TMPSTRUCT.pre_stim{i}.raw(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.All_preceding.([statfield '_fudgeDST'])(end-NumTrials2+1:end);
        TMPSTRUCT.pre_stim{i}.times(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.All_preceding.timevals_fudgeDST_minust0(end-NumTrials2+1:end);
        
        TMPSTRUCT.stim_catch_start{i}.raw(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimCatch.([statfield '_fudgeDST'])(1:NumTrials2);
        TMPSTRUCT.stim_catch_start{i}.times(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimCatch.timevals_fudgeDST_minust0(1:NumTrials2);
        
        TMPSTRUCT.stim_notcatch_start{i}.raw(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimNotCatch.([statfield '_fudgeDST'])(1:NumTrials2);
        TMPSTRUCT.stim_notcatch_start{i}.times(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimNotCatch.timevals_fudgeDST_minust0(1:NumTrials2);
        
        TMPSTRUCT.stim_catch_end{i}.raw(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimCatch.([statfield '_fudgeDST'])(end-NumTrials2+1:end);
        TMPSTRUCT.stim_catch_end{i}.times(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimCatch.timevals_fudgeDST_minust0(end-NumTrials2+1:end);
        
        TMPSTRUCT.stim_notcatch_end{i}.raw(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimNotCatch.([statfield '_fudgeDST'])(end-NumTrials2+1:end);
        TMPSTRUCT.stim_notcatch_end{i}.times(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.StimNotCatch.timevals_fudgeDST_minust0(end-NumTrials2+1:end);
        
        try % might not have enough rends, or not have any data
            TMPSTRUCT.post_stim{i}.raw(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.All_post.([statfield '_fudgeDST'])(1:NumTrials2);
            TMPSTRUCT.post_stim{i}.times(:,k)=StimEpochs_aligned.timewindow{i}.epoch{k}.data.All_post.timevals_fudgeDST_minust0(1:NumTrials2);
        catch err
            TMPSTRUCT.post_stim{i}.raw(:,k)=nan(NumTrials2,1);
            TMPSTRUCT.post_stim{i}.times(:,k)=nan(NumTrials2,1);
        end
        
    end
    
    
    %     collect stats
    TrialFields=fieldnames(TMPSTRUCT);
    for j=1:length(TrialFields);
        trialfield=TrialFields{j};
        
        TMPSTRUCT.(trialfield){i}.mean=mean(TMPSTRUCT.(trialfield){i}.raw);
        TMPSTRUCT.(trialfield){i}.std=std(TMPSTRUCT.(trialfield){i}.raw);
        TMPSTRUCT.(trialfield){i}.sem=std(TMPSTRUCT.(trialfield){i}.raw)/sqrt(NumTrials2-1);
        TMPSTRUCT.(trialfield){i}.meantime=mean(TMPSTRUCT.(trialfield){i}.times);
    end
end


% PLOT - not relative to relative to real time
for i=1:length(TimeFieldsOfInterest);
    lt_figure; hold on;
    
    numtrials=length(TMPSTRUCT.pre_stim{i}.mean);
    
    for ii=1:numtrials
        subplot(ceil(sqrt(numtrials)),ceil(sqrt(numtrials)),ii);
        hold on; title(['Stim epoch ' num2str(ii)]);
        
        errorbar(1, TMPSTRUCT.pre_stim{i}.mean(ii), TMPSTRUCT.pre_stim{i}.sem(ii),'ok');
        
        errorbar(2, TMPSTRUCT.stim_catch_start{i}.mean(ii), TMPSTRUCT.stim_catch_start{i}.sem(ii),'or');
        errorbar(3, TMPSTRUCT.stim_catch_end{i}.mean(ii), TMPSTRUCT.stim_catch_end{i}.sem(ii),'or');
        
        errorbar(2, TMPSTRUCT.stim_notcatch_start{i}.mean(ii), TMPSTRUCT.stim_notcatch_start{i}.sem(ii),'og');
        errorbar(3, TMPSTRUCT.stim_notcatch_end{i}.mean(ii), TMPSTRUCT.stim_notcatch_end{i}.sem(ii),'og');
        
        try % somestimes poststim is missing.
            errorbar(4, TMPSTRUCT.post_stim{i}.mean(ii), TMPSTRUCT.post_stim{i}.sem(ii),'ok');
            % draw a line between pre-stim and post-stim
            line([1 4], [TMPSTRUCT.pre_stim{i}.mean(ii) TMPSTRUCT.post_stim{i}.mean(ii)],'LineStyle','--','Color','b')
        catch err
        end
    end
    lt_subtitle(['N trials: ' num2str(NumTrials2) '; time window ' num2str(i)]);
end

% PLOT - relative to real time
for i=1:length(TimeFieldsOfInterest);
    lt_figure; hold on;
    
    numtrials=length(TMPSTRUCT.pre_stim{i}.mean);
    
    for ii=1:numtrials
        subplot(ceil(sqrt(numtrials)),ceil(sqrt(numtrials)),ii);
        hold on; title(['stim epoch ' num2str(ii)]);
        
        errorbar(TMPSTRUCT.pre_stim{i}.meantime(ii), TMPSTRUCT.pre_stim{i}.mean(ii), TMPSTRUCT.pre_stim{i}.sem(ii),'ok');
        
        errorbar(TMPSTRUCT.stim_catch_start{i}.meantime(ii), TMPSTRUCT.stim_catch_start{i}.mean(ii), TMPSTRUCT.stim_catch_start{i}.sem(ii),'or');
        errorbar(TMPSTRUCT.stim_catch_end{i}.meantime(ii), TMPSTRUCT.stim_catch_end{i}.mean(ii), TMPSTRUCT.stim_catch_end{i}.sem(ii),'or');
        
        errorbar(TMPSTRUCT.stim_notcatch_start{i}.meantime(ii), TMPSTRUCT.stim_notcatch_start{i}.mean(ii), TMPSTRUCT.stim_notcatch_start{i}.sem(ii),'og');
        errorbar(TMPSTRUCT.stim_notcatch_end{i}.meantime(ii), TMPSTRUCT.stim_notcatch_end{i}.mean(ii), TMPSTRUCT.stim_notcatch_end{i}.sem(ii),'og');
        
        try % somestimes poststim is missing.
            errorbar(TMPSTRUCT.post_stim{i}.meantime(ii), TMPSTRUCT.post_stim{i}.mean(ii), TMPSTRUCT.post_stim{i}.sem(ii),'ok');
            % draw a line between pre-stim and post-stim
            line([TMPSTRUCT.pre_stim{i}.meantime(ii) TMPSTRUCT.post_stim{i}.meantime(ii)], [TMPSTRUCT.pre_stim{i}.mean(ii) TMPSTRUCT.post_stim{i}.mean(ii)],'LineStyle','--','Color','b')
        catch err
        end
    end
    lt_subtitle(['Using real time: time window ' num2str(i)]);
end



% CALCULATE DIFFERENCE FROM PRE-STIM AND POST-STIM
% 1) not using temporal information
for i=1:length(TimeFieldsOfInterest);
    numtrials=length(TMPSTRUCT.pre_stim{i}.mean);
    
    TMPSTRUCT2.dev_from_nostim{i}.stim_catch_start=TMPSTRUCT.stim_catch_start{i}.mean-TMPSTRUCT.pre_stim{i}.mean;
    TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_start=TMPSTRUCT.stim_notcatch_start{i}.mean-TMPSTRUCT.pre_stim{i}.mean;
    
    TMPSTRUCT2.dev_from_nostim{i}.stim_catch_end=TMPSTRUCT.post_stim{i}.mean-TMPSTRUCT.stim_catch_end{i}.mean;
    TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_end=TMPSTRUCT.post_stim{i}.mean-TMPSTRUCT.stim_notcatch_end{i}.mean;
    
    % PLOT
    figure;
    % plot without taking absolute value.
    subplot(2,1,1); hold on;
    title(['Deviation from no-stim; time window: ' num2str(i)]);
    
    plot(1,TMPSTRUCT2.dev_from_nostim{i}.stim_catch_start,'or');
    plot(2,TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_start,'og');
    
    plot(4,TMPSTRUCT2.dev_from_nostim{i}.stim_catch_end,'or');
    plot(5,TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_end,'og');
    
    % plot means
    plot(1.2,mean(TMPSTRUCT2.dev_from_nostim{i}.stim_catch_start),'sk','MarkerFaceColor','r');
    plot(2.2,mean(TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_start),'sk','MarkerFaceColor','g');
    plot(4.2,mean(TMPSTRUCT2.dev_from_nostim{i}.stim_catch_end),'sk','MarkerFaceColor','r');
    plot(5.2,mean(TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_end),'sk','MarkerFaceColor','g');
    
    xlim([0.7 5.7]);
    line(xlim, [0 0]);
    
    % plot absolute value
    
    subplot(2,1,2); hold on;
    title(['Absoluted deviation from no-stim; time window: ' num2str(i)]);
    
    plot(1,abs(TMPSTRUCT2.dev_from_nostim{i}.stim_catch_start),'or');
    plot(2,abs(TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_start),'og');
    
    plot(4,abs(TMPSTRUCT2.dev_from_nostim{i}.stim_catch_end),'or');
    plot(5,abs(TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_end),'og');
    
    % plot means
    plot(1.2,mean(abs(TMPSTRUCT2.dev_from_nostim{i}.stim_catch_start)),'sk','MarkerFaceColor','r');
    plot(2.2,mean(abs(TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_start)),'sk','MarkerFaceColor','g');
    plot(4.2,mean(abs(TMPSTRUCT2.dev_from_nostim{i}.stim_catch_end)),'sk','MarkerFaceColor','r');
    plot(5.2,mean(abs(TMPSTRUCT2.dev_from_nostim{i}.stim_notcatch_end)),'sk','MarkerFaceColor','g');
    
    
    xlim([0.7 5.7]);
    
end


% Calculate deviation from line connecting pre and post (i.e. using
% temporal information
tmpfields={'stim_catch_start','stim_notcatch_start','stim_catch_end','stim_notcatch_end'};

for i=1:length(TimeFieldsOfInterest);
    numepochs=length(TMPSTRUCT.pre_stim{i}.mean);
    
    
    post_minus_pre_ff=TMPSTRUCT.post_stim{i}.mean-TMPSTRUCT.pre_stim{i}.mean;
    post_minus_pre_time=TMPSTRUCT.post_stim{i}.meantime-TMPSTRUCT.pre_stim{i}.meantime;
    
    slopes(:,i)=post_minus_pre_ff./post_minus_pre_time;
    
    % find out expected value each each of 4 stim timepoints
    for j=1:length(tmpfields);
        tmpf=tmpfields{j};
        
        TMPSTRUCT.(tmpf){i}.times_from_pre_stim=TMPSTRUCT.(tmpf){i}.meantime-TMPSTRUCT.pre_stim{i}.meantime;
        TMPSTRUCT.(tmpf){i}.expected_mean_if_no_stim_effects=TMPSTRUCT.pre_stim{i}.mean'+TMPSTRUCT.(tmpf){i}.times_from_pre_stim'.*slopes(:,i);
    end
end

% PLOT DEVIATIONS FROM EXPECTED VALUE
% CALCULATE DIFFERENCE FROM PRE-STIM AND POST-STIM
% 1) not using temporal information
for i=1:length(TimeFieldsOfInterest);
    numtrials=length(TMPSTRUCT.pre_stim{i}.mean);
    
    inds=~isnan(TMPSTRUCT.stim_catch_start{i}.expected_mean_if_no_stim_effects); % avoid nans
    TMPSTRUCT2.dev_from_expectval{i}.stim_catch_start=TMPSTRUCT.stim_catch_start{i}.mean(inds)'-TMPSTRUCT.stim_catch_start{i}.expected_mean_if_no_stim_effects(inds);
    TMPSTRUCT2.dev_from_expectval{i}.stim_catch_end=TMPSTRUCT.stim_catch_end{i}.mean(inds)'-TMPSTRUCT.stim_catch_end{i}.expected_mean_if_no_stim_effects(inds);
    TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_start=TMPSTRUCT.stim_notcatch_start{i}.mean(inds)'-TMPSTRUCT.stim_notcatch_start{i}.expected_mean_if_no_stim_effects(inds);
    TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_end=TMPSTRUCT.stim_notcatch_end{i}.mean(inds)'-TMPSTRUCT.stim_notcatch_end{i}.expected_mean_if_no_stim_effects(inds);
    
    
    % PLOT
    figure;
    % plot without taking absolute value.
    subplot(2,1,1); hold on;
    title(['Deviation from expected value; time window: ' num2str(i)]);
    
    plot(1,TMPSTRUCT2.dev_from_expectval{i}.stim_catch_start,'or');
    plot(2,TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_start,'og');
    
    plot(4,TMPSTRUCT2.dev_from_expectval{i}.stim_catch_end,'or');
    plot(5,TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_end,'og');
    
    % plot means
    plot(1.2,mean(TMPSTRUCT2.dev_from_expectval{i}.stim_catch_start),'sk','MarkerFaceColor','r');
    plot(2.2,mean(TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_start),'sk','MarkerFaceColor','g');
    plot(4.2,mean(TMPSTRUCT2.dev_from_expectval{i}.stim_catch_end),'sk','MarkerFaceColor','r');
    plot(5.2,mean(TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_end),'sk','MarkerFaceColor','g');
    
    xlim([0.7 5.7]);
    line(xlim, [0 0]);
    
    % plot absolute value
    
    subplot(2,1,2); hold on;
    title(['Absoluted deviation from expected value; time window: ' num2str(i)]);
    
    plot(1,abs(TMPSTRUCT2.dev_from_expectval{i}.stim_catch_start),'or');
    plot(2,abs(TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_start),'og');
    
    plot(4,abs(TMPSTRUCT2.dev_from_expectval{i}.stim_catch_end),'or');
    plot(5,abs(TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_end),'og');
    
    % plot means
    plot(1.2,mean(abs(TMPSTRUCT2.dev_from_expectval{i}.stim_catch_start)),'sk','MarkerFaceColor','r');
    plot(2.2,mean(abs(TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_start)),'sk','MarkerFaceColor','g');
    plot(4.2,mean(abs(TMPSTRUCT2.dev_from_expectval{i}.stim_catch_end)),'sk','MarkerFaceColor','r');
    plot(5.2,mean(abs(TMPSTRUCT2.dev_from_expectval{i}.stim_notcatch_end)),'sk','MarkerFaceColor','g');
    
    
    xlim([0.7 5.7]);
    
end




