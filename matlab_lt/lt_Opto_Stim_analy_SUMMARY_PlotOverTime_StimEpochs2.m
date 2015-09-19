function [StimEpochs_aligned, Y_day_level_regression]=lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs2(StimEpochs_aligned, PARAMS, DATSTRUCT, NumStimRends, TimeFieldsOfInterest, statfield)
%% PARAMS:

BaselineDays = PARAMS.global.BaselineDays;

% num trials to take to calculate reversion (at start of stim)

% note: will take all data for the day to calculate day stats.


%% REGRESSIONS, DAY LEVELS STATS

RelativeDaysToCompare=[-1]; % i.e. compare to previous day (-1); compare to mean of previous two days ([-2 -1]);
DaysToKeep=1:1000; % days you want to analyze

% Each stim epoch as a data point.
for i=1:length(TimeFieldsOfInterest);
    
    NumEpochs=length(StimEpochs_aligned.timewindow{i}.epoch);
    
    for j=1:NumEpochs;
        
        % ignore baseline days
        if ~any(StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day==BaselineDays);
            
            % ignore any other days you want to ignore
            if any(StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day==DaysToKeep);
                
                % -- FOR EACH EPOCH, COLLECT STATS
                % only look at epochs that have a preceding non-stim window (i.e. to keep
                % things independent
                if isfield(StimEpochs_aligned.timewindow{i}.epoch{j}.data,'All_preceding');
                    
                    % 1) learning (mean pitch of preceding window compared to previous days)
                    % what days to use to determine recent learning?
                    CurDay=StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day;
                    DaysToCompare=RelativeDaysToCompare+CurDay;
                    
                    % 2) get mean of the previous days
                    prevdayvals=[];
                    prevdayvals_time=[];
                    for k=1:length(DaysToCompare);
                        dayind=DaysToCompare(k);
                        
                        % skip days without data
                        %                     if ~isempty(DATSTRUCT_2.data{dayind})
                        %                         tmpvals=DATSTRUCT_2.data{dayind}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST']);
                        %                         tmpvals_time=DATSTRUCT_2.data{dayind}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST;
                        
                        if ~isempty(DATSTRUCT.data{dayind})
                            if isfield(DATSTRUCT.data{dayind}.timewindow{i},'All');
                                tmpvals=DATSTRUCT.data{dayind}.timewindow{i}.All.MINUSBaseHrBins.([statfield '_fudgeDST']);
                                tmpvals_time=DATSTRUCT.data{dayind}.timewindow{i}.All.MINUSBaseHrBins.timevals_fudgeDST;
                                
                                % collect previous days vals
                                prevdayvals=[prevdayvals tmpvals];
                                prevdayvals_time=[prevdayvals_time tmpvals_time];
                            end
                        end
                    end
                    
                    % Collect stats of previous day vals;
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.previous_days.RelativeDaysToCompare=RelativeDaysToCompare;
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.previous_days.DaysToCompare=DaysToCompare;
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.previous_days.([statfield '_fudgeDST'])=prevdayvals;
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.previous_days.timevals_fudgeDST=prevdayvals_time;
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.previous_days.([statfield '_fudgeDST_mean'])=mean(prevdayvals);
                    
                    
                    % 3) get slope of current day (entire day)
                    % values of trials outside stim epoch today
                    if isfield(DATSTRUCT.data{CurDay}.timewindow{i},'All');
                        Yvals=DATSTRUCT.data{CurDay}.timewindow{i}.All.MINUSBaseHrBins.([statfield '_fudgeDST']);
                        Tvals=DATSTRUCT.data{CurDay}.timewindow{i}.All.MINUSBaseHrBins.timevals_fudgeDST;
                        
                        
                        % slope by time
                        X=[ones(length(Tvals),1) Tvals'];
                        Y=Yvals';
                        b=regress(Y,X);
                        
                        StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.current_day.([statfield '_fudgeDST']).slope=b(2);
                        
                        % slope by rendition
                        X=[ones(1,length(Tvals)); 1:length(Tvals)]';
                        Y=Yvals';
                        b=regress(Y,X);
                        
                        StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.current_day.([statfield '_fudgeDST']).slope_rends=b(2);
                        
                    end

                    
                    
                end
            end
        end
    end
end


% -- PERFORM REGRESSIONS
% Collect values into a matrix
clear 'Y_day_level_regression';
for i=1:length(TimeFieldsOfInterest);
    for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
        
        if any(StimEpochs_aligned.timewindow{i}.epoch{j}.epoch_info.day==DaysToKeep)
            if isfield(StimEpochs_aligned.timewindow{i}.epoch{j}.data,'All_preceding');
                
                % reversion (stim notcatch minus stim catch)
                Y_day_level_regression{i}.StimEpochAsDatapoint.reversion(j)=mean(StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.([statfield '_fudgeDST']))-...
                    mean(StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.([statfield '_fudgeDST']));
                
                % recent learning: stim epoch pre-stim mean minus preceding days mean.
                Y_day_level_regression{i}.StimEpochAsDatapoint.recent_learning_rel_previous_days(j)=mean(StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST']))-...
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.previous_days.ffvals_fudgeDST_mean;
                
                % slope of current day (non-stim epochs) (time)
                Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_time(j)=...
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.current_day.([statfield '_fudgeDST']).slope;
                
                % slope of current day (non-stim epochs) (rends)
                Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_rends(j)=...
                    StimEpochs_aligned.timewindow{i}.epoch{j}.stats.day_level_regressions.current_day.([statfield '_fudgeDST']).slope_rends;
                
%                 % Slope of current day (non-stim epoch before stim) (time)
%                 Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_prestim_time(j)= ...
%                     StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope;
%                 
%                  % Slope of current day (non-stim epoch before stim) (rend)
%                 Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_prestim_rend(j)= ...
%                     StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend;
               
            end
        end
    end
end


% == PLOT REGRESSIONS
figure; hold on;
for i=1:length(TimeFieldsOfInterest);
    subplot(2,2,i); hold on;
    % X is recent learning (compared to previous days)
    X=Y_day_level_regression{i}.StimEpochAsDatapoint.recent_learning_rel_previous_days';
    
    % Y is reversion
    Y=Y_day_level_regression{i}.StimEpochAsDatapoint.reversion';
    
    % perform regression
    X=[ones(length(X),1) X];
    [b,~,~,~,stats]=regress(Y,X);
    Slope=b(2);
    
    % plot
    title(['Timefield ' num2str(i) ]);
    xlabel('Recent learning (pre-stim relative to previous days)');
    ylabel('Reversion (hz)');
    plot(X(:,2),Y,'o');
    plot(xlim,b(1) + b(2).*xlim);
    Xlim=xlim;
    Ylim=ylim;
    text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
    
    line(xlim,[0 0],'Color','k','LineStyle','--');
    line([0 0],ylim,'Color','k','LineStyle','--');
    
end
subtitle('Reversion vs. Recent learning (pre-stim minus previous days)')


% slope of current day
figure; hold on;
for i=1:length(TimeFieldsOfInterest);
    subplot(2,2,i); hold on;
    % X is recent learning (compared to previous days)
    X=Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_time';
    
    % Y is reversion
    Y=Y_day_level_regression{i}.StimEpochAsDatapoint.reversion';
    
    % perform regression
    X=[ones(length(X),1) X];
    [b,~,~,~,stats]=regress(Y,X);
    Slope=b(2);
    
    % plot
    title(['Timefield ' num2str(i) ]);
    xlabel('Current day slope (hz/hr)');
    ylabel('Reversion (hz)');
    plot(X(:,2),Y,'o');
    plot(xlim,b(1) + b(2).*xlim);
    Xlim=xlim;
    Ylim=ylim;
    text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
    
    line(xlim,[0 0],'Color','k','LineStyle','--');
    line([0 0],ylim,'Color','k','LineStyle','--');
    
end
subtitle('Reversion vs. Current day slope (hz/hr)')


% % SLOPE OF pre-stim (limited by designated number of rends in bin) (time)
% lt_figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is slope)
%     X=Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_prestim_time';
%     
%     % Y is reversion
%     Y=Y_day_level_regression{i}.StimEpochAsDatapoint.reversion';
%     
%     % perform regression
%     [b,bint,r,rint,stats,SummaryStats]=lt_regress(Y,X,1,0);
%     
%     % plot
%     title(['Timefield ' num2str(i) ]);
%     xlabel('Pre-stim slope (hz/hr)');
%     ylabel('Reversion (hz)');
% %     plot(X(:,2),Y,'o');
% %     plot(xlim,b(1) + b(2).*xlim);
% %     Xlim=xlim;
% %     Ylim=ylim;
% %     text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
% %     
%     line(xlim,[0 0],'Color','k','LineStyle','--');
%     line([0 0],ylim,'Color','k','LineStyle','--');
%     
% end
% lt_subtitle('Reversion vs. pre-stim slope (hz/hr)')
% 
% 
% % SLOPE OF pre-stim (limited by designated number of rends in bin) (rends)
% lt_figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is slope)
%     X=Y_day_level_regression{i}.StimEpochAsDatapoint.current_day_slope_prestim_rend';
%     
%     % Y is reversion
%     Y=Y_day_level_regression{i}.StimEpochAsDatapoint.reversion';
%     
%     % perform regression
%     [b,bint,r,rint,stats,SummaryStats]=lt_regress(Y,X,1,0);
%     
%     % plot
%     title(['Timefield ' num2str(i) ]);
%     xlabel('Pre-stim slope (hz/rend)');
%     ylabel('Reversion (hz)');
% %     plot(X(:,2),Y,'o');
% %     plot(xlim,b(1) + b(2).*xlim);
% %     Xlim=xlim;
% %     Ylim=ylim;
% %     text(Xlim(1)+1,Ylim(1)+1,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))],'Color','r','FontSize',13)
% %     
%     line(xlim,[0 0],'Color','k','LineStyle','--');
%     line([0 0],ylim,'Color','k','LineStyle','--');
%     
% end
% lt_subtitle('Reversion vs. pre-stim slope (hz/rends)')


