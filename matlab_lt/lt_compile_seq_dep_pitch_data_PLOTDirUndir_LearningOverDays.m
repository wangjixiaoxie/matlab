%% PLOT LEARNING OVER DAYS (pitch deviation and zscore)

% UNDIR - Plot figure for this set of fields
figure; hold on;
for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % PLOT MEAN PITCH SUBTRACTING BASELINE
    subplot(length(SylLists.FieldsToPlot),2,-1+j*2); hold on;
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % compile pitch means and SEM for each day, for this syl.
        h=errorbar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
        errorbar_tick(h,200);
    end
    
    legend(FieldsList)
    title(['UNDIR: Mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)')
    xlabel('days')
    
    lt_plot_zeroline
    
    % WN lines
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    
    % PLOT Z-score
    subplot(length(SylLists.FieldsToPlot),2,j*2); hold on;
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % Plot
        plot(DataMatrix.UNDIR.(syl).meanFF_zscore ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
        
    end
    % annotate
    legend(FieldsList)
    title(['UNDIR: Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (Z-scored to baseline)')
    xlabel('days')
    
    lt_plot_zeroline;
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end


% PLOT DIR
if plotDIR==1;
    figure; hold on;
    for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
        FieldsList=SylLists.FieldsToPlot{j};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        % PLOT MEAN PITCH SUBTRACTING BASELINE
        subplot(length(SylLists.FieldsToPlot),2,-1+j*2); hold on;
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % compile pitch means and SEM for each day, for this syl.
            h=errorbar(DataMatrix.DIR.(syl).meanFF_DevFromBase,DataMatrix.DIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
            errorbar_tick(h,200);
        end
        
        legend(FieldsList)
        title(['DIR: Mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (hz) (rel to baseline)')
        xlabel('days')
        
        lt_plot_zeroline
        
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
        
        
        % PLOT Z-score
        subplot(length(SylLists.FieldsToPlot),2,j*2); hold on;
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % Plot
            plot(DataMatrix.DIR.(syl).meanFF_zscore ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
            
        end
        % annotate
        legend(FieldsList)
        title(['DIR: Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (Z-scored to baseline)')
        xlabel('days')
        
        lt_plot_zeroline;
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
end






% PLOT PERCENT DEVIATION FROM BASELINE
% if (0)
%     % PLOT
%     for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
%         FieldsList=SylLists.FieldsToPlot{j};
%
%         plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%
%         % UNDIR - Plot figure for this set of fields
%         figure; hold on;
%         for jj=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%             % Plot
%             plot(DataMatrix.UNDIR.(syl).meanFF_PercentFromBase,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%
%         end
%         % annotate
%         legend(FieldsList)
%         title(['UNDIR: Percent deviation from baseline, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
%         ylabel('FF (% from baseline)')
%         xlabel('days')
%
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%
%
%         % DIR - Plot figure for this set of fields
%         figure; hold on;
%         for jj=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%             % compile pitch means and SEM for each day, for this syl.
%             plot(DataMatrix.DIR.(syl).meanFF_PercentFromBase,'s-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%
%         end
%         %         annotate
%         legend(FieldsList)
%         title(['DIR: Percent deviation from baseline, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
%         ylabel('FF (% from baseline)')
%         xlabel('days')
%
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%     end
%
% end

