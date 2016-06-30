% --- OUTPUTS
Ymeans=[]; % [early late]
Ysems=[];
Yvals={};

% ___________________________________________________________________________________ EARLY DAYS
daystoplot_tmp=Earlydays1;


% ------------------------- Targets
TargetMeans_oneperday=[];
inds_tmp=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds_tmp;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % == Collect each days vals
    TargetMeans_oneperday(daynum)=nanmean(tmp);
end


% -------------- All nontargets
tmp_alldays=[];
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans_oneperday(daynum);
    
    % --- COLLECT
    tmp_alldays=[tmp_alldays tmp']; % columns days, rows=syls
end    

tmp_all_days_vector=reshape(tmp_alldays, numel(tmp_alldays), 1); % put all into one column

% OUTPUT
Yvals{1}=tmp_all_days_vector;
Ymeans(1)=nanmean(tmp_all_days_vector);
Ysems(1)=lt_sem(tmp_all_days_vector);



% ___________________________________________________________________________________ LATE DAYS
daystoplot_tmp=Latedays1;


% ------------------------- Targets
TargetMeans_oneperday=[];
inds_tmp=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds_tmp;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % == Collect each days vals
    TargetMeans_oneperday(daynum)=nanmean(tmp);
end


% -------------- All nontargets
tmp_alldays=[];
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans_oneperday(daynum);
    
    % --- COLLECT
    tmp_alldays=[tmp_alldays tmp']; % columns days, rows=syls
end    

tmp_all_days_vector=reshape(tmp_alldays, numel(tmp_alldays), 1); % put all into one column

% OUTPUT
Yvals{2}=tmp_all_days_vector;
Ymeans(2)=nanmean(tmp_all_days_vector);
Ysems(2)=lt_sem(tmp_all_days_vector);



% _________________________________________________________________________________
% PLOT
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['[All nontargs] Early=' num2str(Earlydays1) '; Late=' num2str(Latedays1)])
xlabel('Early --- Late');
ylabel('fraction of cumulative targ learning')

hbar=lt_plot_bar(1:2, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});
set(hbar, 'FaceColor', 'K');
hold on

p=signrank(Yvals{1}, Yvals{2});
lt_plot_pvalue(p, 'signrank');
