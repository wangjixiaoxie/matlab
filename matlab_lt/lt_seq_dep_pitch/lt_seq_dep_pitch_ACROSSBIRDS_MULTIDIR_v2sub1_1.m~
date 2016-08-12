% ----- one-sided test? (for less than 0) [just for targ 2]
oneSided=1;


% norm to target?
if NormToTarg==1;
    datafield='MeanFFRelBase_NormTargLastBaseDay';
else
    datafield='MeanFFRelBase';
end


if PlotSameTypeOnly==1;
tmp=[DATSTRUCT.information.targ2_sametype_rel_targ1];
ExptInds=find(tmp==1);
elseif PlotSameTypeOnly==2;
    tmp=[DATSTRUCT.information.targ2_sametype_rel_targ1];
ExptInds=find(tmp==0);
elseif PlotSameTypeOnly==3
    % then is same type same seq
    ExptInds=find([DATSTRUCT.information.targ2_sametype_rel_targ1] & [DATSTRUCT.information.targ2_presim_rel_targ1]);
elseif PlotSameTypeOnly==4
    % same type, diff seq
     ExptInds=find([DATSTRUCT.information.targ2_sametype_rel_targ1] & ~[DATSTRUCT.information.targ2_presim_rel_targ1]);
else
    ExptInds=1:length(DATSTRUCT.information);
end


%% === LINE PLOTS + MEAN(SEM)

numexpts=length(DATSTRUCT.data.FirstTarg);
figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];



% ================== targets 1
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('First targets');

Y_all=[];
for i=ExptInds
    Y=DATSTRUCT.data.FirstTarg(i).(datafield);
    X=1:length(Y);
    
    plot(X+0.1, Y, '-k');
    
    % ---- COLLECT
    if ~isempty(Y_all) & size(Y_all,1)~=size(Y, 1)
        disp('SKIPPED AN EXPT as not compatible num of days');
    else
    Y_all=[Y_all Y];
    end
    
end
% --- plot mean
if size(Y_all,2)>1
Ymean=nanmean(Y_all,2);
Ysem=lt_sem(Y_all');
x=1:length(Ymean);


lt_plot(x, Ymean, {'Errors', Ysem, 'LineStyle','-', 'Marker', '+', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8});
end
line(xlim, [0 0], 'LineStyle', '--');


% note sample size:
lt_plot_annotation(1, ['n=' num2str(size(Y_all,2))], 'k')

% ---- line for bidir start
xx=DATSTRUCT.information(1).numPreBidirDays;

line([xx+0.5 xx+0.5], ylim);
% --- significance test for each day (relative to last baseline day)
tmp=Y_all;
day_lastBaseline=DATSTRUCT.information(1).numPreBidirDays;
y_LastBaselineDay=tmp(day_lastBaseline, :);

for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    
    p= signrank(y_thisday, y_LastBaselineDay);
    
    % == plot n and mean(sem)
    nn=numel(y_thisday);
    ymean=mean(y_thisday);
    ysem=lt_sem(y_thisday);
    lt_plot_text(day, 1.2*max(y_thisday), ['n=' num2str(nn) '; mean(sem)=' num2str(ymean) '(' num2str(ysem) ')'], 'm',10);
    
    
    disp([day '-' num2str(p)]);
    if p<0.05;
        lt_plot_text(day, 1.1*max(y_thisday), num2str(p, '%3.2g'), 'r');
    end
    
end
    
% --- significance test for each day (relative to 0)
for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    
    p= signrank(y_thisday);
    
    disp([day '-' num2str(p)]);
    if p<0.15;
        lt_plot_text(day, 1.2*max(y_thisday), [num2str(p, '%3.2g') ' vs 0'], 'g');
    end
    
end
   




% =================== targets 2
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('Second targets');

Y_all=[];
for i=ExptInds
    Y=DATSTRUCT.data.SecondTarg(i).(datafield);
    X=1:length(Y);
    
    plot(X+0.1, Y, '-b');
    
    % ---- COLLECT
        if ~isempty(Y_all) & size(Y_all,1)~=size(Y, 1)
        disp('SKIPPED AN EXPT as not compatible num of days');
    else
    Y_all=[Y_all Y];
    end
end
% --- plot mean
if size(Y_all,2)>1
Ymean=nanmean(Y_all,2);
Ysem=lt_sem(Y_all');
x=1:length(Ymean);


lt_plot(x, Ymean, {'Errors', Ysem, 'LineStyle','-', 'Marker', '+', 'Color', 'b', 'LineWidth', 2, 'MarkerSize', 8});
end
line(xlim, [0 0], 'LineStyle', '--');

% ---- line for bidir start
xx=DATSTRUCT.information(1).numPreBidirDays;

line([xx+0.5 xx+0.5], ylim);
% --- significance test for each day (relative to last baseline day)
tmp=Y_all;
day_lastBaseline=DATSTRUCT.information(1).numPreBidirDays;
y_LastBaselineDay=tmp(day_lastBaseline, :);

for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    
    p= signrank(y_thisday, y_LastBaselineDay);
    
        % == plot n and mean(sem)
    nn=numel(y_thisday);
    ymean=mean(y_thisday);
    ysem=lt_sem(y_thisday);
    lt_plot_text(day, 1.2*max(y_thisday), ['n=' num2str(nn) '; mean(sem)=' num2str(ymean) '(' num2str(ysem) ')'], 'm',10);

    
    disp([day '-' num2str(p)]);
    if p<0.05;
        lt_plot_text(day, 1.1*max(y_thisday), num2str(p, '%3.2g'), 'r');
    end
    
end
% --- significance test for each day (relative to 0)
for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    if oneSided==1
        p= signrank(y_thisday, 0, 'tail', 'left');
        
        disp([day '-' num2str(p)]);
        if p<0.15;
            lt_plot_text(day, 1.2*max(y_thisday), [num2str(p, '%3.2g') ' < 0'], 'g');
        end
        
        
    else
        p= signrank(y_thisday);
        
        disp([day '-' num2str(p)]);
        if p<0.15;
            lt_plot_text(day, 1.2*max(y_thisday), [num2str(p, '%3.2g') ' vs 0'], 'g');
        end
    end
end




%% ===== using boxplots

% ================== targets 1
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('First targets');

Y_all=[];
for i=ExptInds
    Y=DATSTRUCT.data.FirstTarg(i).(datafield);
    X=1:length(Y);
    
%     plot(X+0.1, Y, '-k');
    
    % ---- COLLECT
    if ~isempty(Y_all) & size(Y_all,1)~=size(Y, 1)
        disp('SKIPPED AN EXPT as not compatible num of days');
    else
    Y_all=[Y_all Y];
    end
    
end
% --- plot mean
if size(Y_all,2)>1

boxplot(Y_all', 1:size(Y_all',2), 'color','k');
end
line(xlim, [0 0], 'LineStyle', '--');


% note sample size:
lt_plot_annotation(1, ['n=' num2str(size(Y_all,2))], 'k')

% ---- line for bidir start
xx=DATSTRUCT.information(1).numPreBidirDays;

line([xx+0.5 xx+0.5], ylim);
% --- significance test for each day (relative to last baseline day)
tmp=Y_all;
day_lastBaseline=DATSTRUCT.information(1).numPreBidirDays;
y_LastBaselineDay=tmp(day_lastBaseline, :);

for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    
    p= signrank(y_thisday, y_LastBaselineDay);
    
    % == plot n and mean(sem)
    nn=numel(y_thisday);
    ymean=mean(y_thisday);
    ysem=lt_sem(y_thisday);
    lt_plot_text(day, 1.2*max(y_thisday), ['n=' num2str(nn) '; mean(sem)=' num2str(ymean) '(' num2str(ysem) ')'], 'm',10);
    
    
    disp([day '-' num2str(p)]);
    if p<0.05;
        lt_plot_text(day, 1.1*max(y_thisday), num2str(p, '%3.2g'), 'r');
    end
    
end
    
% --- significance test for each day (relative to 0)
for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    
    p= signrank(y_thisday);
    
    disp([day '-' num2str(p)]);
    if p<0.15;
        lt_plot_text(day, 1.2*max(y_thisday), [num2str(p, '%3.2g') ' vs 0'], 'g');
    end
    
end
   




% =================== targets 2
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('Second targets');

Y_all=[];
for i=ExptInds
    Y=DATSTRUCT.data.SecondTarg(i).(datafield);
    X=1:length(Y);
    
%     plot(X+0.1, Y, '-b');
    
    % ---- COLLECT
        if ~isempty(Y_all) & size(Y_all,1)~=size(Y, 1)
        disp('SKIPPED AN EXPT as not compatible num of days');
    else
    Y_all=[Y_all Y];
    end
end
% --- plot mean
if size(Y_all,2)>1

boxplot(Y_all', 1:size(Y_all',2), 'color','b');
end
line(xlim, [0 0], 'LineStyle', '--');

% ---- line for bidir start
xx=DATSTRUCT.information(1).numPreBidirDays;

line([xx+0.5 xx+0.5], ylim);
% --- significance test for each day (relative to last baseline day)
tmp=Y_all;
day_lastBaseline=DATSTRUCT.information(1).numPreBidirDays;
y_LastBaselineDay=tmp(day_lastBaseline, :);

for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    
    p= signrank(y_thisday, y_LastBaselineDay);
    
        % == plot n and mean(sem)
    nn=numel(y_thisday);
    ymean=mean(y_thisday);
    ysem=lt_sem(y_thisday);
    lt_plot_text(day, 1.2*max(y_thisday), ['n=' num2str(nn) '; mean(sem)=' num2str(ymean) '(' num2str(ysem) ')'], 'm',10);

    
    disp([day '-' num2str(p)]);
    if p<0.05;
        lt_plot_text(day, 1.1*max(y_thisday), num2str(p, '%3.2g'), 'r');
    end
    
end
% --- significance test for each day (relative to 0)
for day=1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase_NormTargLastBaseDay);
    
    y_thisday=tmp(day,:);
    if oneSided==1
        p= signrank(y_thisday, 0, 'tail', 'left');
        
        disp([day '-' num2str(p)]);
        if p<0.15;
            lt_plot_text(day, 1.2*max(y_thisday), [num2str(p, '%3.2g') ' < 0'], 'g');
        end
        
        
    else
        p= signrank(y_thisday);
        
        disp([day '-' num2str(p)]);
        if p<0.15;
            lt_plot_text(day, 1.2*max(y_thisday), [num2str(p, '%3.2g') ' vs 0'], 'g');
        end
    end
end

