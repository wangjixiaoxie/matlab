% 1) For 3-day binned
for i=1:length(SylLists.FieldsInOrder);
    FieldsList=SylLists.FieldsInOrder{i};
    for ii=1:length(FieldsList);
        syl=FieldsList{ii};
        for iii=1:RunWinInds; % how many day bins are there?
            
            % UNDIR
            % 1) learning (hz min baseli)
            EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline(iii,ii)=...
                EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
            
            % 2) generalization
            for j=1:length(SylLists.TargetSyls);
                targsyl=SylLists.TargetSyls{j};
                
                EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                    EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                
            end
            
            
            % DIR
            % 1) learning (hz min baseli)
            EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline(iii,ii)=...
                EpochData.WNdaysSlidingWin{iii}.DIR.(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
            
            % 2) generalization
            for j=1:length(SylLists.TargetSyls);
                targsyl=SylLists.TargetSyls{j};
                
                EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                    EpochData.WNdaysSlidingWin{iii}.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
            end
        end
    end
    EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder=FieldsList;
    EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder=FieldsList;
    
end

% 2) For single-days
for i=1:length(SylLists.FieldsInOrder);
    FieldsList=SylLists.FieldsInOrder{i};
    for ii=1:length(FieldsList);
        syl=FieldsList{ii};
        
        % UNDIR
        % 1) learning (hz min baseli)
        DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline(:,ii)=...
            DataMatrix.UNDIR.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        % 2) generalization
        for j=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{j};
            
            DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
                DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz;
        end
        
        
        % DIR
        % 1) learning (hz min baseli)
        DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline(:,ii)=...
            DataMatrix.DIR.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        % 2) generalization
        for j=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{j};
            
            DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
                DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz;
        end
    end
    DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder=FieldsList;
    DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder=FieldsList;
    
end



% SKIP PLOTTING IN 3D - hard to visualize, rather use heat map or overlayed
% plots

% LEARNING (Hz from baseline) - Then plot in 3d
% for i=1:length(SylLists.FieldsInOrder);
%     figure; hold on;
%
%     % UNDIR
%     subplot(1,2,1); hold on;
%     %     imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
%     %     plot3(ones(11,1),1:11,EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline(:,1));
%     bar3(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
%
%     title('UNDIR');
%     %     colormap('gray'); colorbar
%     %     ylabel(['Binned days during WN, Binsize: ' num2str(RunWind) ', slid 1 day']);
%     ylabel('FF (minus baseline, Hz)');
%     Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%     set(gca,'XTick',Xtick);
%     set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%
%
%     % DIR
%     subplot(1,2,2); hold on;
%     bar3(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
%
%     title('DIR');
%     %     colormap('gray'); colorbar
%     ylabel('FF (minus baseline, Hz)');
%     Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%     set(gca,'XTick',Xtick);
%     set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%     %
%
%     % global title
%     %     subtitle('Learning (mean Hz minus baseline) for each syllable over days.');
% end


% UNDIR
for i=1:length(SylLists.FieldsInOrder);
    figure; hold on;
    
    X=size(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline,2); % how many syls?
    Y=EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline;
    Ymin=min(min(Y));
    Ymax=max(max(Y));
    Ylimits=[Ymin-10 Ymax+10];
    Ylimits_zoom=Ylimits./5;
    
    for ii=1:X;
        subplot(X,2,(-1+(X-ii+1)*2)); hold on;
        bar(Y(:,ii));
        ylim(Ylimits)
        ylabel(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder{ii});
        
        if ii==X;
            xlabel('3-day bin #');
        end
    end
    
    
    for ii=1:X;
        subplot(X,2,((X-ii+1)*2)); hold on;
        bar(Y(:,ii));
        ylim(Ylimits_zoom)
        ylabel(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder{ii});
        if ii==X;
            xlabel('3-day bin #');
        end
        
    end
    subtitle('UNDIR: mean FF in 3-day bin (minus baseline, hz)');
end

% DIR
if plotDIR==1;
    for i=1:length(SylLists.FieldsInOrder);
        figure; hold on;
        
        X=size(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline,2); % how many syls?
        Y=EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline;
        Ymin=min(min(Y));
        Ymax=max(max(Y));
        yy=max(abs([Ymin Ymax]));
        Ylimits=[-yy-10 yy+10];
        Ylimits_zoom=Ylimits./5;
        
        for ii=1:X;
            subplot(X,2,(-1+(X-ii+1)*2)); hold on;
            bar(Y(:,ii));
            ylim(Ylimits)
            ylabel(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder{ii});
            
            if ii==X;
                xlabel('3-day bin #');
            end
        end
        
        
        for ii=1:X;
            subplot(X,2,((X-ii+1)*2)); hold on;
            bar(Y(:,ii));
            ylim(Ylimits_zoom)
            ylabel(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder{ii});
            if ii==X;
                xlabel('3-day bin #');
            end
            
        end
        subtitle('DIR: mean FF in 3-day bin (minus baseline, hz)');
    end
end


% PLOT EACH DAY INDIVIDUALLY

% UNDIR
for i=1:length(SylLists.FieldsInOrder);
    FieldsList=SylLists.FieldsInOrder{i};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    figure; hold on;
    for ii=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{ii}; % actual syl name (e.g. 'a')
        subplot(length(FieldsList),2,-1+(length(FieldsList)-ii+1)*2); hold on;
        
        % Plot pitch means and SEM for each day, for this syl.
        bar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase);
        
        ylim(Ylimits)
        ylabel(syl);
        if ii==length(FieldsList);
            xlabel('3-day bin #');
        end
        
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
    
    % ZOOM IN
    for ii=1:length(FieldsList);
        syl=FieldsList{ii}; % actual syl name (e.g. 'a')
        subplot(length(FieldsList),2,(length(FieldsList)-ii+1)*2); hold on;
        bar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase);
        
        ylim(Ylimits_zoom)
        ylabel(syl);
        if ii==length(FieldsList);
            xlabel('3-day bin #');
        end
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
    
    subtitle('UNDIR - (Syls in order) Mean learning per day (Hz from baseline)');
end


% using zscore
% ylimits=[];
% for i=1:length(SylLists.FieldsInOrder);
%     FieldsList=SylLists.FieldsInOrder{i};
%
%     plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%     figure; hold on;
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         subplot(length(FieldsList),2,-1+ii*2); hold on;
%
%         % Plot pitch means and SEM for each day, for this syl.
%         bar(DataMatrix.UNDIR.(syl).meanFF_zscore);
%
%         ylimits=[ylimits; ylim];
%
%         ylabel(syl);
%         if ii==length(FieldsList);
%             xlabel('3-day bin #');
%         end
%
%         % WN lines
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%     end
%
%     Ylim1=max(max(ylimits));
%     Ylim2=min(max(ylimits));
%     ylimits_zoom=[Ylim2 Ylim1]./5;
%     % ZOOM IN
%     for ii=1:length(FieldsList);
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         subplot(length(FieldsList),2,ii*2); hold on;
%         bar(DataMatrix.UNDIR.(syl).meanFF_zscore);
%
%         ylim(ylimits_zoom)
%         ylabel(syl);
%         if ii==length(FieldsList);
%             xlabel('3-day bin #');
%         end
%         % WN lines
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%     end
%
%     subtitle('UNDIR - (Syls in order) Mean learning per day (Hz from baseline)');
% end




% DIR
if plotDIR==1;
    for i=1:length(SylLists.FieldsInOrder);
        FieldsList=SylLists.FieldsInOrder{i};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        figure; hold on;
        for ii=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{ii}; % actual syl name (e.g. 'a')
            subplot(length(FieldsList),2,-1+(length(FieldsList)-ii+1)*2); hold on;
            
            % Plot pitch means and SEM for each day, for this syl.
            bar(DataMatrix.DIR.(syl).meanFF_DevFromBase);
            
            ylim(Ylimits)
            ylabel(syl);
            if ii==length(FieldsList);
                xlabel('3-day bin #');
            end
            
            % WN lines
            Fn_AnnotateWNLines(plotWNdays,ylim)
        end
        
        % ZOOM IN
        for ii=1:length(FieldsList);
            syl=FieldsList{ii}; % actual syl name (e.g. 'a')
            subplot(length(FieldsList),2,(length(FieldsList)-ii+1)*2); hold on;
            bar(DataMatrix.DIR.(syl).meanFF_DevFromBase);
            
            ylim(Ylimits_zoom)
            ylabel(syl);
            if ii==length(FieldsList);
                xlabel('3-day bin #');
            end
            % WN lines
            Fn_AnnotateWNLines(plotWNdays,ylim)
        end
        
        subtitle('DIR - (Syls in order) Mean learning per day (Hz from baseline)');
    end
end





% OLD METHOD - USING LINE PLOT, NOT BARS.  BARS BETTER.
% PLOT for each day. Put syls in
% order and make subplots
% for i=1:length(SylLists.FieldsInOrder);
%     FieldsList=SylLists.FieldsInOrder{i};
%
%     plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%     figure; hold on;
%
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%
%         subplot(length(FieldsList),2,-1+ii*2); hold on;
%
%         % Plot pitch means and SEM for each day, for this syl.
%         h=errorbar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{ii},'MarkerSize',12);
%         errorbar_tick(h,200);
%
%         ylim(Ylimits)
%         title(syl);
%         ylabel('Pitch (minus baseline) (hz)');
%
%         % WN lines
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%
%         line(xlim,[0 0],'LineStyle','--');
%
%
%     end
%     for ii=1:length(FieldsList);
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         subplot(length(FieldsList),2,ii*2); hold on;
%         h=errorbar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{ii},'MarkerSize',12);
%         errorbar_tick(h,200);
%
%         ylim(Ylimits_zoom)
%         title(syl);
%         ylabel('Pitch (minus baseline) (hz)');
%
%
%         % WN lines
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%
%         line(xlim,[0 0],'LineStyle','--');
%     end
%
%
%     subtitle('UNDIR - Mean Pitch Over Days');
% end


% DO THE SAME FOR DIR
% for i=1:length(SylLists.FieldsInOrder);
%     FieldsList=SylLists.FieldsInOrder{i};
%
%     plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%     figure; hold on;
%
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%
%         subplot(length(FieldsList),2,-1+ii*2); hold on;
%
%         % Plot pitch means and SEM for each day, for this syl.
%         h=errorbar(DataMatrix.DIR.(syl).meanFF_DevFromBase,DataMatrix.DIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{ii},'MarkerSize',12);
%         errorbar_tick(h,200);
%
%         ylim(Ylimits)
%         title(syl);
%         ylabel('Pitch (minus baseline) (hz)');
%
%         % WN lines
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%
%         line(xlim,[0 0],'LineStyle','--');
%
%
%     end
%     for ii=1:length(FieldsList);
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         subplot(length(FieldsList),2,ii*2); hold on;
%         h=errorbar(DataMatrix.DIR.(syl).meanFF_DevFromBase,DataMatrix.DIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{ii},'MarkerSize',12);
%         errorbar_tick(h,200);
%
%         ylim(Ylimits_zoom)
%         title(syl);
%         ylabel('Pitch (minus baseline) (hz)');
%
%
%         % WN lines
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%
%         line(xlim,[0 0],'LineStyle','--');
%     end
%
%
%     subtitle('DIR - Mean Pitch Over Days');
% end






% BELOW - PLOTTING (3D, 3-day bins)

%     title('UNDIR');
%     %     colormap('gray'); colorbar
%     %     ylabel(['Binned days during WN, Binsize: ' num2str(RunWind) ', slid 1 day']);
%     ylabel('FF (minus baseline, Hz)');
%     Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%     set(gca,'XTick',Xtick);
%     set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%
%
%     % DIR
%     subplot(1,2,2); hold on;
%     bar3(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
%
%     title('DIR');
%     %     colormap('gray'); colorbar
%     ylabel('FF (minus baseline, Hz)');
%     Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%     set(gca,'XTick',Xtick);
%     set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%     %
%
%     % global title
%     %     subtitle('Learning (mean Hz minus baseline) for each syllable over days.');
%





% GENERALIZATION (% learning at target) -
% 3-d plot, 3-day bins
% for iii=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{iii};
%
%     for i=1:length(SylLists.FieldsInOrder);
%         figure; hold on;
%
%         % UNDIR
%         subplot(1,2,1); hold on;
%         %     imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl));
%         bar3(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl));
%
%         title('UNDIR');
%         %     colormap('gray'); colorbar
%         %     ylabel(['Binned days during WN, Binsize: ' num2str(RunWind) ', slid 1 day']);
%         ylabel(['Generalization from' targsyl]);
%         Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%         set(gca,'XTick',Xtick);
%         set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%         zlim([-0.5 0.5])
%
%         % DIR
%         subplot(1,2,2); hold on;
%         bar3(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl));
%
%         title('DIR');
%         %     colormap('gray'); colorbar
%         %     ylabel(['Binned days during WN, Binsize: ' num2str(RunWind) ', slid 1 day']);
%         ylabel(['Generalization from' targsyl]);
%         Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%         set(gca,'XTick',Xtick);
%         set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%         zlim([-0.5 0.5])
%
%         % global title
%         %     subtitle('Generalization (% learning rel to target) for each syllable over days.');
%     end
% end
