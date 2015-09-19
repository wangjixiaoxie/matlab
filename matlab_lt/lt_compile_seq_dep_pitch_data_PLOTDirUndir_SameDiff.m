% PLOT LEARNING in hz and zscore.
figure; hold on;
% 1) learning as hz
subplot(1,2,1); hold on;

% PLOT SIMILAR SYLS
FieldsList=SylLists.SylsSame;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii)];
    end
    
    plot(Y,'-ob');
    
    Xlims=xlim;
    text(Xlims(2)+0.1,Y(2),syl,'Color','b'); % label the datapoint
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','b');

% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);




% PLOT DIFF SYLS
FieldsList=SylLists.SylsDifferent;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii)];
    end
    
    plot(Y,'-or');
    Xlims=xlim;
    text(Xlims(2)+0.15,Y(2),syl,'Color','r'); % label the datapoint
    
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','r');

% Format plot
xlim([Xlims(1)-0.2 Xlims(2)+0.4])

ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);
title('FF (minus baseline): Same (RED) and Different (BLUE) syls');
lt_plot_zeroline;



% PLOT SAME, BUT Z-score
subplot(1,2,2); hold on;

% PLOT SIMILAR SYLS
FieldsList=SylLists.SylsSame;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(end,ii)];
    end
    
    plot(Y,'-ob');
    
    Xlims=xlim;
    text(Xlims(2)+0.1,Y(2),syl,'Color','b'); % label the datapoint
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','b');

% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);




% PLOT DIFF SYLS
FieldsList=SylLists.SylsDifferent;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(end,ii)];
    end
    
    plot(Y,'-or');
    Xlims=xlim;
    text(Xlims(2)+0.15,Y(2),syl,'Color','r'); % label the datapoint
    
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','r');

% Format plot
xlim([Xlims(1)-0.2 Xlims(2)+0.4])

ylabel('Pitch (z-score)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);
title('Z-score: Same (RED) and Different (BLUE) syls');
lt_plot_zeroline;


subtitle('Mean values from binned days - first and last days of WN');








% WILL NOT PLOT ABSOLUTE VALUES

% % PLOT AS ABOVE, but absolute values.
% figure; hold on;
% % 1) learning as hz
% subplot(1,2,1); hold on;
% 
% % PLOT SIMILAR SYLS
% FieldsList=SylLists.SylsSame;
% YY=[];
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).meanFF_minusBaseline),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline)];
%     
%     plot(X,Y,'-ob');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
% 
% title('Similar Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% 
% % PLOT DIFF SYLS
% FieldsList=SylLists.SylsDifferent;
% YY=[];
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).meanFF_minusBaseline),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline)];
%     
%     plot(X,Y,'-or');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
% 
% 
% title('Different Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% 
% % zero line
% line(xlim,[0 0],'Color','k','LineStyle','--')
% 
% % global title
% title('ABSOLUTE LEARNING (hz minus baseline) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
% 
% 
% 
% % PLOT SAME, BUT Z-score
% subplot(1,2,2); hold on;
% 
% % PLOT SIMILAR SYLS
% FieldsList=SylLists.SylsSame;
% YY=[];
% 
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).Zscore),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore)];
%     
%     plot(X,Y,'-ob');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
% 
% 
% title('Similar Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% 
% % PLOT DIFF SYLS
% FieldsList=SylLists.SylsDifferent;
% YY=[];
% 
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).Zscore),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore)];
%     
%     plot(X,Y,'-or');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
% 
% 
% title('Different Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% line(xlim,[0 0],'Color','k','LineStyle','--')
% 
% % global title
% title('ABSOLUTE LEARNING (Z-scored) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
% 



% WILL NOT PLOT GENERALIZATION - only meaningful when the target is not
% changing pitch
% PLOT SAME, BUT Generalization
% figure; hold on;
% 
% for i=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{i};
%     subplot(length(SylLists.TargetSyls),1,i); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl),...
%             EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)];
%         
%         plot(X,Y,'-ob');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl),...
%             EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)];
%         
%         plot(X,Y,'-or');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     subtitle(['GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
% end



% PLOT SAME, BUT ABSOLUTE Generalization
% figure; hold on;
% 
% for i=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{i};
%     subplot(length(SylLists.TargetSyls),1,i); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl)),...
%             abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl))];
%         
%         plot(X,Y,'-ob');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl)),...
%             abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl))];
%         
%         plot(X,Y,'-or');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     subtitle(['Absolute GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
% end


%% EXAACTLY LIKE ABOVE, BUT DIR
% NEED TO DO, below is not latest code

% if plotDIR==1;
%     % PLOT LEARNING in hz and zscore.
%     figure; hold on;
%     % 1) learning as hz
%     subplot(1,2,1); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline];
%         
%         plot(X,Y,'-ob');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline];
%         
%         plot(X,Y,'-or');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % zero line
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('LEARNING (hz minus baseline) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     % PLOT SAME, BUT Z-score
%     subplot(1,2,2); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore];
%         
%         plot(X,Y,'-ob');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore];
%         
%         plot(X,Y,'-or');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('LEARNING (Z-scored) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     % PLOT AS ABOVE, but absolute values.
%     figure; hold on;
%     % 1) learning as hz
%     subplot(1,2,1); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline)];
%         
%         plot(X,Y,'-ob');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline)];
%         
%         plot(X,Y,'-or');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % zero line
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('ABSOLUTE LEARNING (hz minus baseline) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     
%     
%     
%     
%     
%     % PLOT SAME, BUT Z-score
%     subplot(1,2,2); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore)];
%         
%         plot(X,Y,'-ob');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore)];
%         
%         plot(X,Y,'-or');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('ABSOLUTE LEARNING (Z-scored) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     
%     
%     % PLOT SAME, BUT Generalization
%     figure; hold on;
%     
%     for i=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{i};
%         subplot(length(SylLists.TargetSyls),1,i); hold on;
%         
%         % PLOT SIMILAR SYLS
%         FieldsList=SylLists.SylsSame;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl),...
%                 EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)];
%             
%             plot(X,Y,'-ob');
%             text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%         
%         
%         title('Similar Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         % PLOT DIFF SYLS
%         FieldsList=SylLists.SylsDifferent;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl),...
%                 EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)];
%             
%             plot(X,Y,'-or');
%             text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%         
%         title('Different Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--')
%         
%         % global title
%         subtitle(['GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
%     end
%     
%     
%     
%     % PLOT SAME, BUT ABSOLUTE Generalization
%     figure; hold on;
%     
%     for i=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{i};
%         subplot(length(SylLists.TargetSyls),1,i); hold on;
%         
%         % PLOT SIMILAR SYLS
%         FieldsList=SylLists.SylsSame;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl)),...
%                 abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl))];
%             
%             plot(X,Y,'-ob');
%             text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%         
%         
%         title('Similar Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         % PLOT DIFF SYLS
%         FieldsList=SylLists.SylsDifferent;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl)),...
%                 abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl))];
%             
%             plot(X,Y,'-or');
%             text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%         
%         title('Different Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--')
%         
%         % global title
%         subtitle(['Absolute GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
%     end
%     
%     
% end
% 

%% SAME AS ABOVE, BUT ALL DAYS, NOT BINNED DAYS


% PLOT LEARNING in hz and zscore.
figure; hold on;


% PLOT SIMILAR SYLS
subplot(1,2,1); hold on;
FieldsList=SylLists.SylsSame;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;

    
    plot(DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(:,ii),'-ob')
    text(NumDays+0.4,DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii),syl,'Color','b'); % label the datapoint
end

% plot mean and SEM
Y=DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_MinusBase.mean;
Ysd=DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_MinusBase.SD;
Ysem=DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_MinusBase.SD/...
    sqrt(length(FieldsList)-1);

errorbar(Y,Ysem,'-ok','MarkerFaceColor','b');

% Plot the target syl in diff color
targsyl=SylLists.TargetSyls{1};
    plot(DataMatrix.UNDIR.(targsyl).meanFF_DevFromBase,'-og')
    text(NumDays+0.4,DataMatrix.UNDIR.(targsyl).meanFF_DevFromBase(end),targsyl,'Color','k'); % label the datapoint



% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
xlabel('days');

Ylimits=ylim;
Fn_AnnotateWNLines(plotWNdays,ylim);
lt_plot_zeroline;



% PLOT DIFF SYLS
subplot(1,2,2); hold on;
FieldsList=SylLists.SylsDifferent;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;

    
    plot(DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(:,ii),'-or')
    text(NumDays+0.4,DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii),syl,'Color','r'); % label the datapoint
end

% plot mean and SEM
Y=DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_MinusBase.mean;
Ysd=DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_MinusBase.SD;
Ysem=DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_MinusBase.SD/...
    sqrt(length(FieldsList)-1);

errorbar(Y,Ysem,'-ok','MarkerFaceColor','r');
ylim(Ylimits)

% Format plot
title('Different Syllables');
ylabel('Pitch (minus baseline) (hz)');
xlabel('days');

Fn_AnnotateWNLines(plotWNdays,ylim);
lt_plot_zeroline;


