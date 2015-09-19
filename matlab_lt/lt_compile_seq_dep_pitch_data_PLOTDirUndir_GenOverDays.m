% PLOT GENERALIZATION SCORE OVER DAYS

% UNDIR - Plot figure for this set of fields
for kk=1:length(SylLists.TargetSyls);
    targsyl=SylLists.TargetSyls{kk};
    
    % Line plot over days.
    % One plot for dir and one for undir
    figure; hold on;
    for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
        FieldsList=SylLists.FieldsToPlot{j};
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        
        subplot(length(SylLists.FieldsToPlot),2,2*j-1); hold on;
        
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            if strcmp(syl,targsyl)==0; % don't plot the target.
                % Plot
                plot(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
            end
        end
        % annotate
        legend(FieldsList)
        title(['UNDIR: Generalization (using Hz) from target : ' targsyl ', ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('Generalization (learning (hz) divided by target syl learning)')
        xlabel('days')
        ylim([-1 1]);
        
        lt_plot_zeroline;
        
        Fn_AnnotateWNLines(plotWNdays,ylim)
        
        
        
        % SAME, but with z-score
        subplot(length(SylLists.FieldsToPlot),2,2*j); hold on;
        
        % UNDIR - Plot figure for this set of fields
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            if strcmp(syl,targsyl)==0; % don't plot the target.
                % Plot
                plot(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
            end
        end
        % annotate
        legend(FieldsList)
        title(['UNDIR: Generalization (Using Z-score) from target : ' targsyl ', ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('Generalization (learning (hz) divided by target syl learning)')
        xlabel('days')
        ylim([-1 1]);
        
        lt_plot_zeroline;
        
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
end






if plotDIR==1;
    % DIR - Plot figure for this set of fields
    for kk=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{kk};
        
        % Line plot over days.
        % One plot for dir and one for undir
        figure; hold on;
        for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
            FieldsList=SylLists.FieldsToPlot{j};
            plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
            
            
            subplot(length(SylLists.FieldsToPlot),2,2*j-1); hold on;
            
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                if strcmp(syl,targsyl)==0; % don't plot the target.
                    % Plot
                    plot(DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
                end
            end
            % annotate
            legend(FieldsList)
            title(['DIR: Generalization (using Hz) from target : ' targsyl ', ' num2str(FirstDay) ' to ' num2str(LastDay)]);
            ylabel('Generalization (learning (hz) divided by target syl learning)')
            xlabel('days')
            ylim([-1 1]);
            
            lt_plot_zeroline;
            
            Fn_AnnotateWNLines(plotWNdays,ylim)
            
            
            
            % SAME, but with z-score
            subplot(length(SylLists.FieldsToPlot),2,2*j); hold on;
            
            % UNDIR - Plot figure for this set of fields
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                if strcmp(syl,targsyl)==0; % don't plot the target.
                    % Plot
                    plot(DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingZ ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
                end
            end
            % annotate
            legend(FieldsList)
            title(['DIR: Generalization (Using Z-score) from target : ' targsyl ', ' num2str(FirstDay) ' to ' num2str(LastDay)]);
            ylabel('Generalization (learning (hz) divided by target syl learning)')
            xlabel('days')
            ylim([-1 1]);
            
            lt_plot_zeroline;
            
            Fn_AnnotateWNLines(plotWNdays,ylim)
        end
    end
end


% if (0)
%     % PLOT generalization scores that are averaged over days for last three days
%     for kk=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{kk};
%         figure; hold on;
%
%         % Get all syls
%         FieldsList=SylFieldsAll;
%         plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%
%
%         % UNDIR - Plot figure for this set of fields
%         subplot(2,1,1); hold on;
%         for jj=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%             % Plot
%             plot(jj,LearningEpochs.UNDIR.(syl).LastThreeDays.GeneralizationFrom.(targsyl).UsingHz.meanFF,...
%                 'o-','Color','k','MarkerFaceColor','k','MarkerSize',12);
%
%
%         end
%
%         % annotate
%         Xtick=1:length(FieldsList); % one tick for each syl. needed.
%         set(gca,'XTick',Xtick);
%         set(gca,'XTickLabel',SylFieldsAll)
%         title(['UNDIR: Generalization from ' targsyl ', average (day = datapoint) over last 3 WN days']);
%         ylabel('Generalization (learning (hz) divided by target syl learning)')
%         ylim([-0.5 0.5]);
%         line(xlim,[0 0],'LineStyle','--')
%
%
%         % DIR - Plot figure for this set of fields
%         subplot(2,1,2); hold on;
%
%         for jj=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%             % compile pitch means and SEM for each day, for this syl.
%             plot(jj,LearningEpochs.DIR.(syl).LastThreeDays.GeneralizationFrom.(targsyl).UsingHz.meanFF,...
%                 'o-','Color','k','MarkerFaceColor','k','MarkerSize',12);
%         end
%
%         %         annotate
%         Xtick=1:length(FieldsList); % one tick for each syl. needed.
%         set(gca,'XTick',Xtick);
%         set(gca,'XTickLabel',SylFieldsAll)
%         title(['DIR: Generalization from ' targsyl ', average (day = datapoint) over last 3 WN days']);
%         ylabel('Generalization (learning (hz) divided by target syl learning)')
%         ylim([-0.5 0.5]);
%         line(xlim,[0 0],'LineStyle','--')
%
%     end
% end

