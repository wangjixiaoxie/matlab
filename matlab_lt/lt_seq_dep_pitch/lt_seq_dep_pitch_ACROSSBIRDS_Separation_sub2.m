% ---- RUN
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['Early=' num2str(Earlydays) '; Late=' num2str(Latedays)])
xlabel('Early --- Late');
ylabel('shift, fraction of target, avg over days')


% --- OUTPUTS
Ymeans=[]; % [early late]
Ysems=[];
Yvals={};

% ___________________________________________________________________________________ EARLY DAYS
daystoplot_tmp=Earlydays;

% -- Targets
TargetMeans_oneperday=[];
inds_targs=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
    
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds_targs;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    TargetMeans_oneperday(daynum)=nanmean(tmp);   
end



% ---  NONTARGETS
tmp_alldays=[];
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
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
daystoplot_tmp=Latedays;

% -- Targets
TargetMeans_oneperday=[];
inds_targs=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
    
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds_targs;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    TargetMeans_oneperday(daynum)=nanmean(tmp);   
end



% ---  NONTARGETS
tmp_alldays=[];
for i=daystoplot_tmp;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
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
% similar
hbar=lt_plot_bar(1:2, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});
set(hbar, 'FaceColor', 'K');
hold on

% sign rank
p=signrank(Yvals{1}, Yvals{2});
lt_plot_pvalue(p, 'signrank');

% % ttest
% [h, p]=ttest(Yvals{1}, Yvals{2});
% lt_plot_pvalue(p, 'signrank');
