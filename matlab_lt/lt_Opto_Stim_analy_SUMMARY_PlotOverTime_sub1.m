

%% LT 3/26/15 - only makes sense in conte4xt of calling program (lt_Opto_Stim_analy_SUMMARY_PlotOverTime)

% THINGS TO DO:
% 1) small windows
% 2) parametrize bin sizes and locations
% 3) pre/post
% 4) WN direction changes
% 5) multilinear regression
% 6) confirm accurate data



% 
% 
% %% FOR ALL STIM EPOCHS, plot stats - reversion versus various things
% 
% for i=1:length(TimeFieldsOfInterest);
%     
%     NumEpochs=length(StimEpochs_aligned.timewindow{i}.epoch);
%     
%     for j=1:NumEpochs;
%         
%         if isfield(StimEpochs_aligned.timewindow{i}.epoch{j}.data,'All_preceding');
%             
%             % COLLECT STATS on this epoch
%             
%             % 1) preceding (nonstim)
%             n=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST'])); % num renditions
%             if ~isempty(NumTrials_otherstats);
%                 % if there are not enough data, give use warning
%                 if n<NumTrials_otherstats;
%                     disp(['warning, not enough trials for preceding (nonstim) for day: ' num2str(ii), ' epoch: ' num2str(k) ' have ' num2str(length(Inds)) ' trials. Will take all available rends']);
%                     
%                     Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST']);
%                     Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.timevals_fudgeDST_minust0;
%                     
%                 else % have enough trials
%                     
%                     Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST'])(end-NumTrials_otherstats+1:end);
%                     Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.timevals_fudgeDST_minust0((end-NumTrials_otherstats+1:end));
%                     
%                 end
%                 
%             end
%             
%             StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Yvals);
%             
%             
%             %                     StimEpochs_2.data.timewindow{i}.All_preceding.mean(c)=mean(Yvals);
%             
%             % slope
%             % by time
%             X=[ones(length(Tvals),1) Tvals'];
%             Y=Yvals';
%             b=regress(Y,X);
%             
%             %                                         figure; hold on;
%             %                                         plot(X(:,2),Y,'o');
%             %                                         plot(xlim,b(1) + b(2).*xlim);
%             
%             StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope=b(2);
%             
%             % slope by rendition
%             X=[ones(1,length(Tvals)); 1:length(Tvals)]';
%             Y=Yvals';
%             b=regress(Y,X);
%             
%             StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend=b(2);
%             
%             
%             
%             % 2) Stim epoch (catch and not catch)
%             tmpfields={'StimCatch','StimNotCatch'};
%             
%             for kkk=1:length(tmpfields);
%                 tmpf=tmpfields{kkk};
%                 n=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])); % num renditions
%                 if ~isempty(NumTrials_reversion);
%                     if n<NumTrials_reversion;
%                         disp(['warning, not enough trials for ' tmpf ' for day: ' num2str(ii), ' epoch: ' num2str(k) ' have ' num2str(length(Inds)) ' trials. Will continue']);
%                         
%                         Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST']);
%                         Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0;
%                         
%                     else % have enough trials
%                         
%                         Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])(1:NumTrials_reversion);
%                         Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0(1:NumTrials_reversion);
%                     end
%                 end
%                 
%                 
%                 StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Yvals);
%                 
%                 % slope
%                 % by time
%                 X=[ones(length(Tvals),1) Tvals'];
%                 Y=Yvals';
%                 b=regress(Y,X);
%                 
%                 %                     figure; hold on;
%                 %                     plot(X(:,2),Y,'o');
%                 %                     plot(xlim,b(1) + b(2).*xlim);
%                 
%                 StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset.([statfield '_fudgeDST']).slope=b(2);
%                 
%                 % by rendition
%                 X=[ones(1,length(Tvals)); 1:length(Tvals)]';
%                 Y=Yvals';
%                 b=regress(Y,X);
%                 
%                 StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend=b(2);
%             end
%         end
%     end
% end
% 
% 
% 
% %% PLOT REGRESSION (stim effects vs. recent pitch change
% 
% % first extract data into format: cell array, each cell containing one time
% % field, and that contains matrix where rows are trials and columns go:
% 
% % for comparing mean values
% for i=1:length(TimeFieldsOfInterest);
%     for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
%         Ymean{i}(j,1)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
%         Ymean{i}(j,2)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
%         Ymean{i}(j,3)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
%     end
% end
% 
% % for slopes
% for i=1:length(TimeFieldsOfInterest);
%     for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
%         Yslope{i}(j,1)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope;
%         Yslope{i}(j,2)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
%         Yslope{i}(j,3)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
%         Yslope{i}(j,4)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend;
%         Yslope{i}(j,5)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend;
%         Yslope{i}(j,6)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset.([statfield '_fudgeDST']).slope;
%     end
% end
% 
% 
% % PLOT
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
%     RealStats.timewindow{i}.regressions.recent_learning.p=stats(3);
%     RealStats.timewindow{i}.regressions.recent_learning.r2=stats(1);
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. Recent learning (StimCatch - pre-stim)')
% 
% % PLOT REGRESSION (stim effects vs. current level of learning (try for both
% % catch and pre-stim)
% figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is mean of pre-stim
%     X=Ymean{i}(:,1);
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
%     title(['Timefield '  num2str(i) ]);
%     xlabel('Learning (prestim)');
%     ylabel('Reversion (StimNotCatch minus StimCatch)');
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
%     RealStats.timewindow{i}.regressions.current_learning_prestim.p=stats(3);
%     RealStats.timewindow{i}.regressions.current_learning_prestim.r2=stats(1);
%     
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. learning (Pre-stim)')
% 
% 
% % PLOT REGRESSION (stim effects vs. current level of learning (try for both
% % catch and pre-stim)
% figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is mean of stim catch
%     X=Ymean{i}(:,2);
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
%     title(['Timefield '  num2str(i) ]);
%     xlabel('Learning (Stim Catch)');
%     ylabel('Reversion (StimNotCatch minus StimCatch)');
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
%     RealStats.timewindow{i}.regressions.current_learning_stimcatch.p=stats(3);
%     RealStats.timewindow{i}.regressions.current_learning_stimcatch.r2=stats(1);
%     
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. learning (StimCatch)')
% 
% 
% 
% % PLOT
% figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is slope of pre-stim
%     X=Yslope{i}(:,1);
%     
%     % Y is reversion (stim not catch minus stim catch)
%     Y=diff(Yslope{i},1,2);
%     Y=Y(:,2);
%     
%     % perform regression
%     X=[ones(length(X),1) X];
%     [b,~,~,~,stats]=regress(Y,X);
%     Slope=b(2);
%     
%     % plot
%     title(['Timefield '  num2str(i) ]);
%     xlabel('Slope of pre-stim (hz/hr)');
%     ylabel('Reversion (StimNotCatch minus StimCatch)');
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
%     RealStats.timewindow{i}.regressions.slope_time_prestim.p=stats(3);
%     RealStats.timewindow{i}.regressions.slope_time_prestim.r2=stats(1);
%     
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (time) of pre-stim')
% 
% 
% 
% % PLOT REGRESSION (during stim slope (stim catch) (by time)) -----------------------
% % extract data to matrix
% 
% % PLOT
% figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is slope of stim catch
%     X=Yslope{i}(:,6);
%     
%     % Y is reversion (stim not catch minus stim catch)
%     Y=diff(Yslope{i},1,2);
%     Y=Y(:,2);
%     
%     % perform regression
%     X=[ones(length(X),1) X];
%     [b,~,~,~,stats]=regress(Y,X);
%     Slope=b(2);
%     
%     % plot
%     title(['Timefield '  num2str(i) ]);
%     xlabel('Slope of Durign stim (catch) (hz/hr)');
%     ylabel('Reversion (StimNotCatch minus StimCatch)');
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
%     RealStats.timewindow{i}.regressions.slope_time_stimcatch.p=stats(3);
%     RealStats.timewindow{i}.regressions.slope_time_stimcatch.r2=stats(1);
%     
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (time) during stim (catch)')
% 
% 
% % PLOT
% figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is slope of pre-stim
%     X=Yslope{i}(:,4);
%     
%     % Y is reversion (stim not catch minus stim catch)
%     Y=diff(Yslope{i},1,2);
%     Y=Y(:,2);
%     
%     % perform regression
%     X=[ones(length(X),1) X];
%     [b,~,~,~,stats]=regress(Y,X);
%     Slope=b(2);
%     
%     % plot
%     title(['Timefield '  num2str(i) ]);
%     xlabel('Slope of pre-stim (hz/rend)');
%     ylabel('Reversion (StimNotCatch minus StimCatch)');
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
%     RealStats.timewindow{i}.regressions.slope_rend_prestim.p=stats(3);
%     RealStats.timewindow{i}.regressions.slope_rend_prestim.r2=stats(1);
%     
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (rends) of pre-stim')
% 
% 
% 
% % PLOT REGRESSION (during stim slope (stim catch) (by rends)) -----------------------
% % PLOT
% figure; hold on;
% for i=1:length(TimeFieldsOfInterest);
%     subplot(2,2,i); hold on;
%     % X is slope of stim catch
%     X=Yslope{i}(:,5);
%     
%     % Y is reversion (stim not catch minus stim catch)
%     Y=diff(Yslope{i},1,2);
%     Y=Y(:,2);
%     
%     % perform regression
%     X=[ones(length(X),1) X];
%     [b,~,~,~,stats]=regress(Y,X);
%     Slope=b(2);
%     
%     % plot
%     title(['Timefield '  num2str(i) ]);
%     xlabel('Slope of Durign stim (catch) (hz/rend)');
%     ylabel('Reversion (StimNotCatch minus StimCatch)');
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
%     RealStats.timewindow{i}.regressions.slope_rend_stimcatch.p=stats(3);
%     RealStats.timewindow{i}.regressions.slope_rend_stimcatch.r2=stats(1);
%     
% end
% subtitle('Reversion (StimNotCatch - StimCatch) vs. Slope (using rends) during stim (catch)')


%% REDO ABOVE, BUT DOING PERMUTATION TESTS
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



