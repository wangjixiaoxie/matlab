%% PLOT MEAN OF LAST FEW DAYS IN ORDER OF SYLLABLES

% PLOT LEARNING AS MINUS BASELINE (HZ);
% UNDIR - Plot figure for this set of fields


Y={};
Ysem={};
for ll=1:length(SylLists.FieldsInOrder);
    % Get all syls
    FieldsList=SylLists.FieldsInOrder{ll};
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    Y{ll}=[];
    Ysem{ll}=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        Y{ll}=[Y{ll} EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline];
        Ysem{ll}=[Ysem{ll} EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).semFF];
    end
    
    Ymax(ll)=max(Y{ll});
    Ymin(ll)=min(Y{ll});
    
end


% Plot
figure; hold on;

for ll=1:length(SylLists.FieldsInOrder);
    subplot(length(SylLists.FieldsInOrder),1,ll); hold on;
    FieldsList=SylLists.FieldsInOrder{ll};
    
    bar(Y{ll})
    
    ylim([min(Ymin)-20 max(Ymax)+20]);
    
    % annotate
    Xtick=1:length(FieldsList); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',FieldsList)
    ylabel('FF (hz), minus baseline')
    
end
subtitle(['UNDIR: Learning (Hz minus baseline) last 3 WN days']);





%         % DIR - Plot figure for this set of fields
%         if plotDIR==1;
%         figure; hold on;
%         for ll=1:length(SylLists.FieldsInOrder);
%             subplot(length(SylLists.FieldsInOrder),1,ll); hold on;
%             % Get all syls
%             FieldsList=SylLists.FieldsInOrder{ll};
%             plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%
%
%             Y=[];
%             for jj=1:length(FieldsList); % how many fields within this set?
%                 syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%                 Y=[Y EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline];
%
%             end
%
%             Ymax=max(Y);
%             Ymin=min(Y);
%
%             % Plot
%             bar(Y);
%
%
%
%             % annotate
%             Xtick=1:length(FieldsList); % one tick for each syl. needed.
%             set(gca,'XTick',Xtick);
%             set(gca,'XTickLabel',FieldsList)
%             ylabel('FF (hz), minus baseline')
%
%         end
%         subtitle(['DIR: Learning (Hz minus baseline) last 3 WN days']);
%         end
%


%% PLOT, but for specified day ranges

% FIRST EXTRACT AND AVERAGE OVER THE DAYS DESIRED

if exist('DaysForSnapshot','var')
    for dd=1:length(DaysForSnapshot);
        days=lt_convert_EventTimes_to_RelTimes(FirstDay,DaysForSnapshot{dd});
        DaysToAverage=days.JustDays_rel(1):days.JustDays_rel(2); % indices
        fnamedays=(['from' DaysForSnapshot{dd}{1} 'to' DaysForSnapshot{dd}{2}]); % field name for structure below
        
        % NOW AVERAGE over those days
        for i=1:length(SylFieldsAll);
            syl=SylFieldsAll{i};
            
            X=[];
            for ii=DaysToAverage;
                try % in case day lacks data
                    X=[X cell2mat(AllDays_compiled_DirAndUndir.UNDIR{ii}.data.(syl)(:,1))']; % collect all data points from baseline days
                catch err
                end
            end
            
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).rawFF=X;
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF=mean(X);
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).stdFF=std(X);
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).n=length(X);
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).semFF=std(X)/(length(X)-1);
            
            
            % 3) DEVIATION FROM BASELINE
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF_minusBaseline=...
                EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF-EpochData.Baseline.UNDIR.(syl).meanFF;
            
            %z-score
            EpochData.Snapshot.(fnamedays).UNDIR.(syl).Zscore=...
                EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF_minusBaseline/EpochData.Baseline.UNDIR.(syl).stdFF;
            
        end
        
        % GENERALIZATION STATS
        for i=1:length(SylFieldsAll);
            syl=SylFieldsAll{i};
            
            for k=1:length(SylLists.TargetSyls);
                targsyl=SylLists.TargetSyls{k}; % target;
                
                EpochData.Snapshot.(fnamedays).UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
                    EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF_minusBaseline/EpochData.Snapshot.(fnamedays).UNDIR.(targsyl).meanFF_minusBaseline;
                
                EpochData.Snapshot.(fnamedays).UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
                    EpochData.Snapshot.(fnamedays).UNDIR.(syl).Zscore/EpochData.Snapshot.(fnamedays).UNDIR.(targsyl).Zscore;
            end
        end
        
        
        % NOW PLOT
        % First, gather data in plotting format for BAR
        Y={};
        Ysem={};
        for ll=1:length(SylLists.FieldsInOrder);
            % Get all syls
            FieldsList=SylLists.FieldsInOrder{ll};
            plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
            
            Y{ll}=[];
            Ysem{ll}=[];
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                Y{ll}=[Y{ll} EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF_minusBaseline];
                Ysem{ll}=[Ysem{ll} EpochData.Snapshot.(fnamedays).UNDIR.(syl).semFF];
            end
            
            Ymax(ll)=max(Y{ll});
            Ymin(ll)=min(Y{ll});
            
        end
        
        
        % Second, Plot
        figure; hold on;
        
        for ll=1:length(SylLists.FieldsInOrder);
            subplot(length(SylLists.FieldsInOrder),1,ll); hold on;
            FieldsList=SylLists.FieldsInOrder{ll};
            
            bar(Y{ll})
            
            ylim([min(Ymin)-20 max(Ymax)+20]);
            
            % annotate
            Xtick=1:length(FieldsList); % one tick for each syl. needed.
            set(gca,'XTick',Xtick);
            set(gca,'XTickLabel',FieldsList)
            ylabel('FF (hz), minus baseline')
            
        end
        subtitle(['UNDIR: Learning (Hz minus baseline) from ' DaysForSnapshot{dd}{1} ' to ' DaysForSnapshot{dd}{2}]);
        
        
    end
end


% NOW PLOT LEARNING


%% PLOT LEARNING AS GENERALIZATION COEFFICIENT - IGNORE, just scales
% learning.  not useful.
%     for kk=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{kk};
%
%         for ll=1:length(SylLists.FieldsInOrder);
%             figure; hold on;
%
%             % Get all syls
%             FieldsList=SylLists.FieldsInOrder{ll};
%             plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%
%
%             % UNDIR - Plot figure for this set of fields
%             subplot(2,1,1); hold on;
%             for jj=1:length(FieldsList); % how many fields within this set?
%                 syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%                 % Plot
%                 plot(jj,EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl),...
%                     'o-','Color','k','MarkerFaceColor','k','MarkerSize',10);
%
%
%             end
%             line(xlim,[0 0],'LineStyle','--')
%
%             % annotate
%             Xtick=1:length(FieldsList); % one tick for each syl. needed.
%             set(gca,'XTick',Xtick);
%             set(gca,'XTickLabel',FieldsList)
%             title(['UNDIR: Generalization from ' targsyl ', average (day = datapoint) over last 3 WN days']);
%             ylabel('Generalization (learning (hz) divided by target syl learning)')
%             ylim([-0.5 0.5]);
%
%
%             % DIR - Plot figure for this set of fields
%             subplot(2,1,2); hold on;
%
%             for jj=1:length(FieldsList); % how many fields within this set?
%                 syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%
%                 plot(jj,EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl),...
%                     's-','Color','k','MarkerFaceColor','k','MarkerSize',12);
%             end
%
%             line(xlim,[0 0],'LineStyle','--')
%
%             %         annotate
%             Xtick=1:length(FieldsList); % one tick for each syl. needed.
%             set(gca,'XTick',Xtick);
%             set(gca,'XTickLabel',FieldsList)
%             title(['DIR: Generalization from ' targsyl ', average (day = datapoint) over last 3 WN days']);
%             ylabel('Generalization (learning (hz) divided by target syl learning)')
%             ylim([-0.5 0.5]);
%         end
%     end

