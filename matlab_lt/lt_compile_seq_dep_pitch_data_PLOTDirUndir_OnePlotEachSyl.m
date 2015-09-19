% PLOT DEVIATION OF MEAN FROM BASELINE
plot_colors=lt_make_plot_colors(length(SylFieldsAll),0); % initiate colors for plot
for j=1:length(SylFieldsAll); % how many fields (i.e. syls)?
    syl=SylFieldsAll{j};
    
    figure; hold on;
    
    % 1) PLOT individual renditions overlayed with mean PITCH (absolute)
    subplot(3,1,1); hold on;
    
    for i=1:NumDays;
        
        % UNDIR
        try % in case no data.
            % get times
            time=cell2mat(AllDays_compiled_DirAndUndir.UNDIR{i}.data.(syl)(:,6));
            [~, B]=lt_convert_datenum_to_hour(time);
            time=B.days+-0.5+i; % so that noon of 1st day is on x=1.
            % get FF
            FF=cell2mat(AllDays_compiled_DirAndUndir.UNDIR{i}.data.(syl)(:,1));
            % Plot
            plot(time,FF,'.','Color','k');
        catch err
        end
        
        % DIR
        try % in case no data.
            % get times
            time=cell2mat(AllDays_compiled_DirAndUndir.DIR{i}.data.(syl)(:,6));
            [~, B]=lt_convert_datenum_to_hour(time);
            time=B.days+-0.5+i; % so that noon of 1st day is on x=1.
            % get FF
            FF=cell2mat(AllDays_compiled_DirAndUndir.DIR{i}.data.(syl)(:,1));
            % Plot
            plot(time,FF,'.','Color','b');
        catch err
        end
        
        
    end
    % Overlay with mean values.
    h=errorbar(DataMatrix.UNDIR.(syl).meanFF,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',8);
    errorbar_tick(h,200);
    
    h=errorbar(DataMatrix.DIR.(syl).meanFF,DataMatrix.UNDIR.(syl).semFF,'s--','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',8);
    errorbar_tick(h,200);
    
    
    %annotate
    legend(SylFieldsAll{j});
    title(['UNDIR (circle) and DIR (square): Mean Pitch, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz)')
    xlabel('days')
    xlim([1 NumDays]);
    
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    % 2) PLOT MEAN PITCH (rel to baseline)
    subplot(3,1,2); hold on;
    
    h=errorbar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase ,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',12);
    errorbar_tick(h,200);
    
    h=errorbar(DataMatrix.DIR.(syl).meanFF_DevFromBase ,DataMatrix.UNDIR.(syl).semFF,'s--','Color',plot_colors{j},'MarkerSize',12);
    errorbar_tick(h,200);
    
    
    %annotate
    %     legend(SylFieldsAll{j});
    title(['UNDIR (circle) and DIR (square): Mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)')
    xlabel('days')
    xlim([1 NumDays]);
    
    % Line for 0
    line(xlim, [0 0],'Color','k','LineStyle','--')
    
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    
    % 3) % PLOT DIFFERENCE BETWEEN DIR AND UNDIR
    subplot(3,1,3); hold on;
    DirMinusUndir=DataMatrix.DIRvsUNDIR.(syl).meanFF_DIRminusUNDIR;
    BaselineDirMinusUndir=EpochData.Baseline.DIRvsUNDIR.(syl).meanFF_DIRminusUNDIR;
    
    PutAFPComp=-(DirMinusUndir-BaselineDirMinusUndir); % i.e. UNDIR minus DIR, taking into account DIR positive pitch effect at baseline.
    
    plot(PutAFPComp,'o-','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',12)
    xlim([1 NumDays]);
    
    % Line for 0
    line(xlim, [0 0],'Color','k','LineStyle','--')
    
    %annotate
    %     legend(SylFieldsAll{j});
    title(['PUTATIVE AFP COMPONENT (UNDIR minus DIR, baseline subtracted): ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('hz')
    xlabel('days')
    
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end

