
% ===== extract data
Y1=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.Edge_rawFF_previous;
Y2=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.Edge_rawFF_current;

X= floor(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.StartTime_SecondEpoch); % days

Xtime_trans=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.StartTime_SecondEpoch; % transition time

% ===== make sure there is only one datapoint per day
if any(diff(X)==0);
    % then a day has more than 1 data
    disp 'MAJOR ERROR more than one transition of a specific type in a day - problem -  - plotting, but check to see what is problem'
end

% ==== make sure all rends in these bins are in the same day (i.e. not
% wraparaound from other day)
Points_contained_within_day=[];
tmp=cell2mat(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ng_Ind1,ng_Ind2}.EpochInds);
epochinds1=tmp(1:2:length(tmp));
epochinds2=tmp(2:2:length(tmp));

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


% ===== get median times for each epoch
tval1_median=[];
tval2_median=[];
for i=1:length(epochinds1);
    epoch1=epochinds1(i);
    epoch2=epochinds2(i);
    
    tval1_median(i)=median(SORTED_DATA.ByNoteGroup(ng_Ind1).Stats_OneDataPtPerEpoch.END.Tvals{epoch1});
    tval2_median(i)=median(SORTED_DATA.ByNoteGroup(ng_Ind2).Stats_OneDataPtPerEpoch.BEGINNING.Tvals{epoch2});
end



% ===== plot each day
for i=1:length(X);
%     day=X(i);
%     
%     if Points_contained_within_day(i)==1;
%         lt_plot([day-0.2 day], [mean(Y1{i}) mean(Y2{i})], {'LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
%     else
%         lt_plot([day-0.2 day], [mean(Y1{i}) mean(Y2{i})], {'MarkerFaceColor','','LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
%     end

if plot_epoch_median_times==1;
    
     % plot at switch times
    xtimes=[tval1_median(i) tval2_median(i)];
    
    % MEAN
%     if Points_contained_within_day(i)==1;
%         lt_plot(xtimes, [mean(Y1{i}) mean(Y2{i})], {'LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
%     else
%         lt_plot(xtimes, [mean(Y1{i}) mean(Y2{i})], {'MarkerFaceColor','none','LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
%     end
   
    if Points_contained_within_day(i)==1;
        lt_plot(xtimes, [median(Y1{i}) median(Y2{i})], {'LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
    else
        lt_plot(xtimes, [median(Y1{i}) median(Y2{i})], {'MarkerFaceColor','none','LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
    end

else
    
    % plot at switch times
    xtime=Xtime_trans(i);
    
    % MEAN
%     if Points_contained_within_day(i)==1;
%         lt_plot([xtime-0.1 xtime+0.1], [mean(Y1{i}) mean(Y2{i})], {'LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
%     else
%         lt_plot([xtime-0.1 xtime+0.1], [mean(Y1{i}) mean(Y2{i})], {'MarkerFaceColor','none','LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
%     end
    
    if Points_contained_within_day(i)==1;
        lt_plot([xtime-0.1 xtime+0.1], [median(Y1{i}) median(Y2{i})], {'LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
    else
        lt_plot([xtime-0.1 xtime+0.1], [median(Y1{i}) median(Y2{i})], {'MarkerFaceColor','none','LineStyle', '-', 'Color', plotcol, 'Errors', [lt_sem(Y1{i}) lt_sem(Y2{i})]});
    end

end
end


