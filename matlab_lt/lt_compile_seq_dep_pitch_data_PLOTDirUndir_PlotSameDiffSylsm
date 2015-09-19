% LEARNING (Hz from baseline) - 3 day bins

% UNDIR
figure; hold on;
for i=1:length(SylLists.FieldsInOrder);
    subplot(1,length(SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Binned days during WN, Binsize: ' num2str(RunWind)]);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
    
    
    %     global title
    subtitle('UNDIR: Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
end

% ZOOM IN
figure; hold on;
for i=1:length(SylLists.FieldsInOrder);
    subplot(1,length(SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Binned days during WN, Binsize: ' num2str(RunWind)]);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
    
    caxis(Ylimits_zoom)
    
    %     global title
    subtitle('UNDIR: Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
end


if plotDIR==1;
    figure; hold on;
    for i=1:length(SylLists.FieldsInOrder);
        subplot(1,length(SylLists.FieldsInOrder),i); hold on;
        
        % UNDIR
        imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
        
        colormap('gray');
        hbar=colorbar;
        ylabel(['Binned days during WN, Binsize: ' num2str(RunWind)]);
        Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
        set(gca,'XTick',Xtick);
        set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
        
        
        %     global title
        subtitle('DIR: Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
    end
    
    % ZOOM IN
    figure; hold on;
    for i=1:length(SylLists.FieldsInOrder);
        subplot(1,length(SylLists.FieldsInOrder),i); hold on;
        
        % UNDIR
        imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
        
        colormap('gray');
        hbar=colorbar;
        ylabel(['Binned days during WN, Binsize: ' num2str(RunWind)]);
        Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
        set(gca,'XTick',Xtick);
        set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
        
        caxis(Ylimits_zoom)
        
        %     global title
        subtitle('DIR: Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
    end
    
end


% PLOT SAME, but day by day
figure; hold on;
for i=1:length(SylLists.FieldsInOrder);
    subplot(1,length(SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Days']);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
    
    % WN start and end line
    line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
    line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
    
    %     global title
    subtitle('UNDIR: Learning (mean Hz minus baseline) for each syllable.');
end

% ZOOM IN

figure; hold on;
for i=1:length(SylLists.FieldsInOrder);
    subplot(1,length(SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Days']);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
    
    caxis(Ylimits_zoom)
    
    % WN start and end line
    
    line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
    line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
    
    %     global title
    subtitle('UNDIR: Learning (mean Hz minus baseline) for each syllable.');
end



% DIR
if plotDIR==1;
    figure; hold on;
    for i=1:length(SylLists.FieldsInOrder);
        subplot(1,length(SylLists.FieldsInOrder),i); hold on;
        
        % UNDIR
        imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
        
        colormap('gray');
        hbar=colorbar;
        ylabel(['Days']);
        Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
        set(gca,'XTick',Xtick);
        set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
        
        % WN start and end line
        line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
        line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
        
        %     global title
        subtitle('UNDIR: Learning (mean Hz minus baseline) for each syllable.');
    end
    
    % ZOOM IN
    
    figure; hold on;
    for i=1:length(SylLists.FieldsInOrder);
        subplot(1,length(SylLists.FieldsInOrder),i); hold on;
        
        % UNDIR
        imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.meanFF_minusBaseline);
        
        colormap('gray');
        hbar=colorbar;
        ylabel(['Days']);
        Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
        set(gca,'XTick',Xtick);
        set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
        
        caxis(Ylimits_zoom)
        
        % WN start and end line
        
        line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
        line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
        
        %     global title
        subtitle('UNDIR: Learning (mean Hz minus baseline) for each syllable.');
    end
end





% % GENERALIZATION (% learning at target) -
% % Then plot in heat map (rows = days), (columns = syl, in order).
% for iii=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{iii};
%
%     for i=1:length(SylLists.FieldsInOrder);
%         figure; hold on;
%
%         % UNDIR
%         subplot(1,2,1); hold on;
%         imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl));
%
%         title('UNDIR');
%             colorbar
%             ylabel(['Binned days during WN, Binsize: ' num2str(RunWind) ', slid 1 day']);
%         ylabel(['Generalization from' targsyl]);
%         Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder); % one tick for each syl. needed.
%         set(gca,'XTick',Xtick);
%         set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder.UNDIR.SylLists.FieldsInOrderNum{i}.FieldNamesInOrder)
%         zlim([-0.5 0.5])
%
%         % DIR
%         subplot(1,2,2); hold on;
%         imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder.DIR.SylLists.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl));
%
%         title('DIR');
%             colorbar
%             ylabel(['Binned days during WN, Binsize: ' num2str(RunWind) ', slid 1 day']);
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
%
