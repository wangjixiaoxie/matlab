function lt_context_PLOT_sub2(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum)

    % ===== extract data
    Y1=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.Edge_rawFF_previous;
    Y2=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.Edge_rawFF_current;
    
    X= floor(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.StartTime_SecondEpoch); % days
%     Xtime_trans=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.StartTime_SecondEpoch; % transition time
    
    tmp=cell2mat(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.EpochInds);
    epochinds1=tmp(1:2:length(tmp));
    epochinds2=tmp(2:2:length(tmp));
    

    % ==== FILTER to keep data for desired phase only
    PhaseNumsEdges=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1, ng_Ind2}.PhaseNum_BothEpochs_EdgeRends;
    inds=PhaseNumsEdges==PhaseNum;
    
%     if isempty(inds)
%         % then no data for this phase
%         break
%     end
    
    Y1=Y1(inds);
    Y2=Y2(inds);
    X=X(inds);
    epochinds1=epochinds1(inds);
    epochinds2=epochinds2(inds);
%     Xtime_trans=Xtime_trans(inds);


    
    % ===== make sure there is only one datapoint per day
    if any(diff(X)==0);
        % then a day has more than 1 data
        disp 'MAJOR ERROR more than one transition of a specific type in a day - problem -  - plotting, but check to see what is problem'
    end

    
    % ==== make sure all rends in these bins are in the same day (i.e. not
    % wraparaound from other day)
    Points_contained_within_day=[];

    for i=1:length(epochinds1);
        epoch1=epochinds1(i);
        epoch2=epochinds2(i);
        
        tvals_endof1=SORTED_DATA.ByNoteGroup(ng_Ind1).Stats_OneDataPtPerEpoch.END.Tvals{epoch1};
        tvals_startof2=SORTED_DATA.ByNoteGroup(ng_Ind2).Stats_OneDataPtPerEpoch.BEGINNING.Tvals{epoch2};
        
        Points_contained_within_day(i)=0; % prime with failure
        % make sure all on same day
        if any(diff(floor(tvals_endof1))) || any(diff(floor(tvals_startof2)))
            disp 'ERROR One epoch does not have enough data to take edge that is strictly within one day - plotting, but check to see what is problem'
        elseif floor(tvals_endof1(1)) ~= floor(tvals_startof2(1));
            disp 'ERROR One transition is exactly overnight - plotting, but check to see what is problem'
        else
            % then is good
            Points_contained_within_day(i)=1;
        end
    end
    
    
%     % ===== get median times for each epoch
%     tval1_median=[];
%     tval2_median=[];
%     for i=1:length(epochinds1);
%         epoch1=epochinds1(i);
%         epoch2=epochinds2(i);
%         
%         tval1_median(i)=median(SORTED_DATA.ByNoteGroup(ng_Ind1).Stats_OneDataPtPerEpoch.END.Tvals{epoch1});
%         tval2_median(i)=median(SORTED_DATA.ByNoteGroup(ng_Ind2).Stats_OneDataPtPerEpoch.BEGINNING.Tvals{epoch2});
%     end
%     
    
    
    % ====== Plot all transitions
    Yall=[];
    for i=1:length(X);
        plot([1 2], [median(Y1{i}) median(Y2{i})], 'o-', 'Color', plotcol);
        
        Yall=[Yall; [median(Y1{i}) median(Y2{i})]];
    end
    
    xlim([0 3]);

    % ===== PLOT MEAN
    if size(Yall,1)>1;
    Ymean=mean(Yall);
    Ysem=lt_sem(Yall);
    lt_plot([1.1 2.1], Ymean, {'Color','k' , 'Errors', Ysem, 'LineStyle', '-'})
    end
    
    
    