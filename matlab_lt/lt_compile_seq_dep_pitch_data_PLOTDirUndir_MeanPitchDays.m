%% PLOT MEAN PITCH OVER DAYS

figure; hold on;

for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % UNDIR - Plot figure for this set of fields
    subplot(length(SylLists.FieldsToPlot),1,j); hold on;
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % Plot
        h=errorbar(DataMatrix.UNDIR.(syl).meanFF,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
        errorbar_tick(h,200);
        
    end
    % annotate
    legend(FieldsList)
    title(['UNDIR: Mean Pitch, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz)')
    xlabel('days')
    
    % WN lines
    Fn_AnnotateWNLines(plotWNdays,ylim)
end


if plotDIR==1;
    figure; hold on;
    for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
        FieldsList=SylLists.FieldsToPlot{j};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        % DIR - Plot figure for this set of fields
        subplot(length(SylLists.FieldsToPlot),1,j); hold on;
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % compile pitch means and SEM for each day, for this syl.
            h=errorbar(DataMatrix.DIR.(syl).meanFF,DataMatrix.DIR.(syl).semFF,'s-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
            errorbar_tick(h,200);
        end
        
        % annotate
        legend(FieldsList)
        title(['DIR: Mean Pitch, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (hz)')
        xlabel('days')
        
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
end
