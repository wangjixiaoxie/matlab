function [StimEpochs_aligned, RealStats] = lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs3(StimEpochs_aligned, StimEpochs, DATSTRUCT, NumTrials_prestim, NumTrials_stim, TimeFieldsOfInterest, statfield, Conditions)

%% PARAMS

    
    % NumTrials_prestim % for stats prestim
    % NumTrials_stim % during stim  (20 means 40 total)
    % RealStats holds regression p and R2;
    % Only_Plot_Days_With_One_Stim = 1, then skips a stim epoch if it is on a day with multiple stims.
    
    
    %% FOR ALL STIM EPOCHS, plot stats - reversion versus various things
    
    for i=1:length(TimeFieldsOfInterest);
        
        NumEpochs=length(StimEpochs_aligned.timewindow{i}.epoch);
        
        for j=1:NumEpochs;
            
            % ==== TEST FILTERING CONDITIONS FOR THIS EPOCH (does it pass?)
            if Conditions.Only_Plot_Days_With_One_Stim==1;
                % check to make sure today only has one stim.
                today_ind=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                num_epochs_today=length(StimEpochs{today_ind}.epoch);
                
                % skip this epoch if today has multiple epochs.
                if num_epochs_today>1;
                    continue
                end
            end
            
            if Conditions.OnlyPlotFirstEpochOfDay==1;
                epochind=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.withinday_epoch;
                if epochind>1;
                    continue;
                end
            end
            
            
            if Conditions.OnlyPlotIfNoStimYesterday==1;
                today_ind=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                
                if isempty(StimEpochs{today_ind-1});
                    DidStimYesterday=0;
                else
                    DidStimYesterday=1;
                end
                
                if DidStimYesterday==1;
                    continue
                end
            end
            % =====================================
            
            
            if isfield(StimEpochs_aligned.timewindow{i}.epoch{j}.data,'All_preceding');
                
                % COLLECT STATS on this epoch
                
                % 1) === preceding (nonstim)
                n=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST'])); % num renditions
                if ~isempty(NumTrials_prestim); % then want to take subset of trials.
                    % if there are not enough data, give user warning
                    if n<NumTrials_prestim;
                        day=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                        withinday_epoch=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.withinday_epoch;
                        disp(['warning, not enough trials for preceding (nonstim) for day: ' num2str(day), ' epoch: ' num2str(withinday_epoch) ' have ' num2str(n) ' trials. Will take all available rends']);
                        
                        Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST']);
                        Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.timevals_fudgeDST_minust0;
                        
                    else % have enough trials
                        
                        Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST'])(end-NumTrials_prestim+1:end);
                        Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.timevals_fudgeDST_minust0((end-NumTrials_prestim+1:end));
                        
                    end
                    
                end
                
                % -- mean
                StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Yvals);
                
                
                % -- slope
                % by time
                X=[ones(length(Tvals),1) Tvals'];
                Y=Yvals';
                b=regress(Y,X);
                
                StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope=b(2);
                
                % slope by rendition
                X=[ones(1,length(Tvals)); 1:length(Tvals)]';
                Y=Yvals';
                b=regress(Y,X);
                
                StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend=b(2);
                
                
                % == MEAN OF ALL OF YESTERDAY
                % get all data for yesterday (average of stim and notstim)
                today_ind=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                
                Dat_Yesterday=[];
                if ~isempty(DATSTRUCT.data{today_ind-1}); % yesterday has data;
                    fname_list=fieldnames(DATSTRUCT.data{today_ind-1}.timewindow{i});
                    
                    for kk=1:length(fname_list);
                        fname=fname_list{kk};
                        
                        Dat_Yesterday=[Dat_Yesterday DATSTRUCT.data{today_ind-1}.timewindow{i}.(fname).MINUSBaseHrBins.([statfield '_fudgeDST'])];
                    end
                end
                
                % -- get yesterday mean ff
                if ~isempty(Dat_Yesterday);
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_yesterday.STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Dat_Yesterday);
                end
                
                
                % 2) === Stim epoch (catch and not catch)
                tmpfields={'StimCatch','StimNotCatch'};
                
                for kkk=1:length(tmpfields);
                    tmpf=tmpfields{kkk};
                    n=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])); % num renditions
                    if ~isempty(NumTrials_stim);
                        if n<NumTrials_stim;
                            day=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                            withinday_epoch=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.withinday_epoch;
                            disp(['warning, not enough trials for ' tmpf ' for day: ' num2str(day), ' epoch: ' num2str(withinday_epoch) ' have ' num2str(n) ' trials. Will take all available rends']);
                            
                            Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST']);
                            Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0;
                            
                        else % have enough trials
                            
                            Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])(1:NumTrials_stim);
                            Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0(1:NumTrials_stim);
                        end
                    end
                    
                    % -- Mean
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Yvals);
                    
                    % -- slope
                    % by time
                    X=[ones(length(Tvals),1) Tvals'];
                    Y=Yvals';
                    b=regress(Y,X);
                    
                    %                     figure; hold on;
                    %                     plot(X(:,2),Y,'o');
                    %                     plot(xlim,b(1) + b(2).*xlim);
                    
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset.([statfield '_fudgeDST']).slope=b(2);
                    
                    % by rendition
                    X=[ones(1,length(Tvals)); 1:length(Tvals)]';
                    Y=Yvals';
                    b=regress(Y,X);
                    
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend=b(2);
                    
                end
                
                % 3) === STIM EPOCH (combining stim and notstim)
                % if asked for N trials, takes first N of stim and N of nonstim
                % trials and combines
                
                tmpfields={'StimCatch','StimNotCatch'};
                Tvals_tot=[];
                Yvals_tot=[];
                
                for kkk=1:length(tmpfields);
                    tmpf=tmpfields{kkk};
                    n=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])); % num renditions
                    if ~isempty(NumTrials_stim);
                        if n<NumTrials_stim;
                            day=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                            withinday_epoch=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.withinday_epoch;
                            disp(['warning, not enough trials for ' tmpf ' for day: ' num2str(day), ' epoch: ' num2str(withinday_epoch) ' have ' num2str(n) ' trials. Will take all available rends']);
                            
                            Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST']);
                            Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0;
                            
                        else % have enough trials
                            
                            Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])(1:NumTrials_stim);
                            Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0(1:NumTrials_stim);
                        end
                    end
                    
                    % == combine stim and notstim data
                    Tvals_tot=[Tvals_tot Tvals];
                    Yvals_tot=[Yvals_tot Yvals];
                    
                end
                
                
                % -- sort by time
                [Tvals_tot, inds] = sort(Tvals_tot);
                Yvals_tot = Yvals_tot(inds);
                
                
                % -- MEAN
                StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimAndNotStim.STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Yvals);
                
            end
        end
    end
    
    
    
    %% PLOT REGRESSION (stim effects vs. recent pitch change
    
    % first extract data into format: cell array, each cell containing one time
    % field, and that contains matrix where rows are trials and columns go:
    
    % for comparing mean values
    for i=1:length(TimeFieldsOfInterest);
        for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
            try
                Ymean{i}(j,1)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
                Ymean{i}(j,2)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
                Ymean{i}(j,3)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
                Ymean{i}(j,4)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimAndNotStim.STATS_PossiblySubset.([statfield '_fudgeDST']).mean; % Stim + Nonstim values, mean.
                try
                    Ymean{i}(j,5)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_yesterday.STATS_PossiblySubset.([statfield '_fudgeDST']).mean; % yesterday non-stim mean
                catch err
                    Ymean{i}(j,5)=nan;
                end
            catch err
            end
        end
    end
    
    % for slopes
    for i=1:length(TimeFieldsOfInterest);
        for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
            try
                Yslope{i}(j,1)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope;
                Yslope{i}(j,2)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
                Yslope{i}(j,3)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
                Yslope{i}(j,4)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend;
                Yslope{i}(j,5)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend;
                Yslope{i}(j,6)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).slope;
            catch err
            end
        end
    end
    
    
    %% =================== PLOT
    
    % == Learning from yesterday (today pre-stim minus yesterday whole day mean)
    lt_figure; hold on;
    
    for i=1:length(TimeFieldsOfInterest);
        lt_subplot(2,2,i); hold on;
        
        X=Ymean{i}(:,1)-Ymean{i}(:,5); % prestim minus yesterday
        
        % Y is Reversion
        Y=diff(Ymean{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        [~,~,~,~,~,SummaryStats] = lt_regress(Y, X, 1);
        
        % plot
        title(['Timefield ' num2str(i) ]);
        xlabel('Learning (today pre-stim minus yesterday mean)');
        ylabel('Reversion (StimNotCatch minus StimCatch)');
        
        line(xlim,[0 0],'Color','k','LineStyle','--');
        line([0 0],ylim,'Color','k','LineStyle','--');
        
        % save stats
        RealStats.timewindow{i}.regressions.learning_today_minus_yesterday.p=SummaryStats.p;
        RealStats.timewindow{i}.regressions.learning_today_minus_yesterday.r2=SummaryStats.R2;
        RealStats.timewindow{i}.regressions.learning_today_minus_yesterday.slope=SummaryStats.slope;
       
    end
    lt_subtitle('Reversion (StimNotCatch - StimCatch) vs. Learning (today pre-stim minus yesterday mean)')
    
    
    
    
    
    % == Recent learning (using stim catch/notcatch mean - preceding)
    lt_figure; hold on;
    
    for i=1:length(TimeFieldsOfInterest);
        lt_subplot(2,2,i); hold on;
        % X is Stimcombined minus preceding
        X=Ymean{i}(:,4)-Ymean{i}(:,1);
        
        % Y is Reversion
        Y=diff(Ymean{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        [~,~,~,~,~,SummaryStats] = lt_regress(Y, X, 1);
        %
        %     X=[ones(length(X),1) X];
        %     [b,~,~,~,stats]=regress(Y,X);
        %     Slope=b(2);
        
        % plot
        title(['Timefield ' num2str(i) ]);
        xlabel('Recent learning (using mean of stim and notstim minus nostim(preceding))');
        ylabel('Reversion (StimNotCatch minus StimCatch)');
        %     lt_plot(X(:,2),Y);
        %     plot(xlim,b(1) + b(2).*xlim);
        %     Xlim=xlim;
        %     Ylim=ylim;
        %     text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
        
        line(xlim,[0 0],'Color','k','LineStyle','--');
        line([0 0],ylim,'Color','k','LineStyle','--');
        
        % save stats
        RealStats.timewindow{i}.regressions.recent_learning_StimCombined.p=SummaryStats.p;
        RealStats.timewindow{i}.regressions.recent_learning_StimCombined.r2=SummaryStats.R2;
        RealStats.timewindow{i}.regressions.recent_learning_StimCombined.slope=SummaryStats.slope;
        
    end
    lt_subtitle('Reversion (StimNotCatch - StimCatch) vs. Recent learning (StimCombined - pre-stim)')
    
    
    % % == Recent learning (using stim catch - preceding)
    % figure; hold on;
    %
    % for i=1:length(TimeFieldsOfInterest);
    %     subplot(2,2,i); hold on;
    %     % X is stim catch minus preceding
    %     X=diff(Ymean{i},1,2);
    %     X=X(:,1);
    %
    %     % Y is stim not catch minus stim catch
    %     Y=diff(Ymean{i},1,2);
    %     Y=Y(:,2);
    %
    %     % perform regression
    %     X=[ones(length(X),1) X];
    %     [b,~,~,~,stats]=regress(Y,X);
    %     Slope=b(2);
    %
    %     % plot
    %     title(['Timefield ' num2str(i) ]);
    %     xlabel('StimCatch minus NoStim(preceding)');
    %     ylabel('StimNotCatch minus StimCatch');
    %     plot(X(:,2),Y,'o');
    %     plot(xlim,b(1) + b(2).*xlim);
    %     Xlim=xlim;
    %     Ylim=ylim;
    %     text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
    %
    %     line(xlim,[0 0],'Color','k','LineStyle','--');
    %     line([0 0],ylim,'Color','k','LineStyle','--');
    %
    %     % save stats
    %     RealStats.timewindow{i}.regressions.recent_learning_StimCatch.p=stats(3);
    %     RealStats.timewindow{i}.regressions.recent_learning_StimCatch.r2=stats(1);
    % end
    % subtitle('Reversion (StimNotCatch - StimCatch) vs. Recent learning (StimCatch - pre-stim)')
    
    
    % PLOT REGRESSION (stim effects vs. current level of learning (try for both
    % catch and pre-stim)
    figure; hold on;
    for i=1:length(TimeFieldsOfInterest);
        subplot(2,2,i); hold on;
        % X is mean of pre-stim
        X=Ymean{i}(:,1);
        
        % Y is stim not catch minus stim catch
        Y=diff(Ymean{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        [~,~,~,~,~,SummaryStats] = lt_regress(Y, X, 1);
        
        % plot
        title(['Timefield '  num2str(i) ]);
        xlabel('Learning (prestim)');
        ylabel('Reversion (StimNotCatch minus StimCatch)');

        % save stats
        RealStats.timewindow{i}.regressions.current_learning_prestim.p=SummaryStats.p;
        RealStats.timewindow{i}.regressions.current_learning_prestim.r2=SummaryStats.R2;
         RealStats.timewindow{i}.regressions.current_learning_prestim.slope=SummaryStats.slope;
       
        
    end
    subtitle('Reversion (StimNotCatch - StimCatch) vs. learning (Pre-stim)')
    
    
%     % PLOT REGRESSION (stim effects vs. current level of learning (try for both
%     % catch and pre-stim)
%     figure; hold on;
%     for i=1:length(TimeFieldsOfInterest);
%         subplot(2,2,i); hold on;
%         % X is mean of stim catch
%         X=Ymean{i}(:,2);
%         
%         % Y is stim not catch minus stim catch
%         Y=diff(Ymean{i},1,2);
%         Y=Y(:,2);
%         
%         % perform regression
%         X=[ones(length(X),1) X];
%         [b,~,~,~,stats]=regress(Y,X);
%         Slope=b(2);
%         
%         % plot
%         title(['Timefield '  num2str(i) ]);
%         xlabel('Learning (Stim Catch)');
%         ylabel('Reversion (StimNotCatch minus StimCatch)');
%         plot(X(:,2),Y,'o');
%         plot(xlim,b(1) + b(2).*xlim);
%         Xlim=xlim;
%         Ylim=ylim;
%         text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--');
%         line([0 0],ylim,'Color','k','LineStyle','--');
%         
%         % save stats
%         RealStats.timewindow{i}.regressions.current_learning_stimcatch.p=stats(3);
%         RealStats.timewindow{i}.regressions.current_learning_stimcatch.r2=stats(1);
%         
%     end
%     subtitle('Reversion (StimNotCatch - StimCatch) vs. learning (StimCatch)')
%     
    
    
    % PLOT
    figure; hold on;
    for i=1:length(TimeFieldsOfInterest);
        subplot(2,2,i); hold on;
        % X is slope of pre-stim
        X=Yslope{i}(:,1);
        
        % Y is reversion (stim not catch minus stim catch)
        Y=diff(Yslope{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        [~,~,~,~,~,SummaryStats] = lt_regress(Y, X, 1);
        
        % plot
        title(['Timefield '  num2str(i) ]);
        xlabel('Slope of pre-stim (hz/hr)');
        ylabel('Reversion (StimNotCatch minus StimCatch)');

        % save stats
        RealStats.timewindow{i}.regressions.slope_time_prestim.p=SummaryStats.p;
        RealStats.timewindow{i}.regressions.slope_time_prestim.r2=SummaryStats.R2;
        RealStats.timewindow{i}.regressions.slope_time_prestim.slope=SummaryStats.slope;
        
    end
    subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (time) of pre-stim')
    
    
    
%     % PLOT REGRESSION (during stim slope (stim catch) (by time)) -----------------------
%     % extract data to matrix
%     
%     % PLOT
%     figure; hold on;
%     for i=1:length(TimeFieldsOfInterest);
%         subplot(2,2,i); hold on;
%         % X is slope of stim catch
%         X=Yslope{i}(:,6);
%         
%         % Y is reversion (stim not catch minus stim catch)
%         Y=diff(Yslope{i},1,2);
%         Y=Y(:,2);
%         
%         % perform regression
%         X=[ones(length(X),1) X];
%         [b,~,~,~,stats]=regress(Y,X);
%         Slope=b(2);
%         
%         % plot
%         title(['Timefield '  num2str(i) ]);
%         xlabel('Slope of Durign stim (catch) (hz/hr)');
%         ylabel('Reversion (StimNotCatch minus StimCatch)');
%         plot(X(:,2),Y,'o');
%         plot(xlim,b(1) + b(2).*xlim);
%         Xlim=xlim;
%         Ylim=ylim;
%         text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--');
%         line([0 0],ylim,'Color','k','LineStyle','--');
%         
%         % save stats
%         RealStats.timewindow{i}.regressions.slope_time_stimcatch.p=stats(3);
%         RealStats.timewindow{i}.regressions.slope_time_stimcatch.r2=stats(1);
%         
%     end
%     subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (time) during stim (catch)')
    
    
    % PLOT
    figure; hold on;
    for i=1:length(TimeFieldsOfInterest);
        subplot(2,2,i); hold on;
        % X is slope of pre-stim
        X=Yslope{i}(:,4);
        
        % Y is reversion (stim not catch minus stim catch)
        Y=diff(Yslope{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        [~,~,~,~,~,SummaryStats] = lt_regress(Y, X, 1);
        
        % plot
        title(['Timefield '  num2str(i) ]);
        xlabel('Slope of pre-stim (hz/rend)');
        ylabel('Reversion (StimNotCatch minus StimCatch)');

        % save stats
        RealStats.timewindow{i}.regressions.slope_rend_prestim.p=SummaryStats.p;
        RealStats.timewindow{i}.regressions.slope_rend_prestim.r2=SummaryStats.R2;
        RealStats.timewindow{i}.regressions.slope_rend_prestim.slope=SummaryStats.slope;
        
    end
    subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (rends) of pre-stim')
    
    
    
%     % PLOT REGRESSION (during stim slope (stim catch) (by rends)) -----------------------
%     % PLOT
%     figure; hold on;
%     for i=1:length(TimeFieldsOfInterest);
%         subplot(2,2,i); hold on;
%         % X is slope of stim catch
%         X=Yslope{i}(:,5);
%         
%         % Y is reversion (stim not catch minus stim catch)
%         Y=diff(Yslope{i},1,2);
%         Y=Y(:,2);
%         
%         % perform regression
%         X=[ones(length(X),1) X];
%         [b,~,~,~,stats]=regress(Y,X);
%         Slope=b(2);
%         
%         % plot
%         title(['Timefield '  num2str(i) ]);
%         xlabel('Slope of Durign stim (catch) (hz/rend)');
%         ylabel('Reversion (StimNotCatch minus StimCatch)');
%         plot(X(:,2),Y,'o');
%         plot(xlim,b(1) + b(2).*xlim);
%         Xlim=xlim;
%         Ylim=ylim;
%         text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--');
%         line([0 0],ylim,'Color','k','LineStyle','--');
%         
%         % save stats
%         RealStats.timewindow{i}.regressions.slope_rend_stimcatch.p=stats(3);
%         RealStats.timewindow{i}.regressions.slope_rend_stimcatch.r2=stats(1);
%         
%     end
%     subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (using rends) during stim (catch)')
    
    
    
    %% === REDO ABOVE, BUT DOING PERMUTATION TESTS
    if (0)
        lt_Opto_Stim_analy_SUMMARY_PlotOverTime_sub1_1
        
        
        % PLOT THE SHUFFLE TEST OUTCOMES
        for i=TimeFieldsOfInterest;
            things_regressed=fieldnames(ShuffleStats.timewindow{i}.regressions);
            
            % Plot histogram of P values
            figure; hold on;
            for j=1:length(things_regressed);
                regressed=things_regressed{j};
                
                subplot(ceil(length(things_regressed)/3),3,j); hold on;
                title(regressed);
                ylabel('counts');
                xlabel('log10(p-value)')
                
                [N, bin]=hist(log10(ShuffleStats.timewindow{i}.regressions.(regressed).p), 100);
                bar(bin,N,'hist');
                
                Ylim=ylim;
                
                % put a line indication alpha values
                Pertiles=prctile(log10(ShuffleStats.timewindow{i}.regressions.(regressed).p),[0.5 5]);
                line([Pertiles(1) Pertiles(1)], ylim);
                line([Pertiles(2) Pertiles(2)], ylim);
                text(Pertiles(1), Ylim(2), 'blue: 0.5 and 5 percentiles')
                
                %         % put a line indication p=0.05
                %         line([log10(0.05) log10(0.05)], ylim,'Color','g');
                %         text(tmp,Ylim(2),
                
                % put a line indication value of real data
                tmp=log10(RealStats.timewindow{i}.regressions.(regressed).p);
                line([tmp tmp], ylim, 'Color','r');
                
                % what is the probability of that real value ocurring in the
                % shuffled data?
                prob_occurance=sum(ShuffleStats.timewindow{i}.regressions.(regressed).p<RealStats.timewindow{i}.regressions.(regressed).p)/...
                    length(ShuffleStats.timewindow{i}.regressions.(regressed).p);
                text(tmp,Ylim(2),['real value, probability: ' num2str(prob_occurance)]);
            end
            subtitle(['Timefield: ' num2str(i) '. Histogram of p-value for regression vs. reversion'])
            
            % Plot histogram of R2
            figure; hold on;
            for j=1:length(things_regressed);
                regressed=things_regressed{j};
                
                subplot(ceil(length(things_regressed)/3),3,j); hold on;
                title(regressed);
                ylabel('counts');
                xlabel('R2')
                
                [N, bin]=hist(ShuffleStats.timewindow{i}.regressions.(regressed).r2,100);
                bar(bin,N,'hist');
                
                % put a line indication alpha values
                Pertiles=prctile(ShuffleStats.timewindow{i}.regressions.(regressed).r2,[95 99.5]);
                line([Pertiles(1) Pertiles(1)], ylim);
                line([Pertiles(2) Pertiles(2)], ylim);
                text(Pertiles(1), Ylim(2), 'blue: 95 and 99.5 percentiles')
                
                % put a line indication value of real data
                tmp=RealStats.timewindow{i}.regressions.(regressed).r2;
                line([tmp tmp], ylim, 'Color','r');
                
                % what is the probability of that real value ocurring in the
                % shuffled data?
                prob_occurance=sum(ShuffleStats.timewindow{i}.regressions.(regressed).r2>RealStats.timewindow{i}.regressions.(regressed).r2)/...
                    length(ShuffleStats.timewindow{i}.regressions.(regressed).r2);
                Ylim=ylim;
                text(tmp,Ylim(2),['real value, probability: ' num2str(prob_occurance)]);
                
            end
            subtitle(['Timefield: ' num2str(i) '. Histogram of R2 for regression vs. reversion'])
        end
        
    end
    
    
