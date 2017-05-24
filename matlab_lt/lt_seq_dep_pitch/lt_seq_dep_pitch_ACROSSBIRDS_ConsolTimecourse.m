function lt_seq_dep_pitch_ACROSSBIRDS_ConsolTimecourse(OUTPUT_any, DATSTRUCT_any)
%% lt 5/8/17 - plotting consolidation experiments smoothed timecourse


%% 


% ===== MULTIDIR
NumBinsToPlot = max(find(~cellfun('isempty', {OUTPUT_any.firsttarget.FFrelBaseInBin.PBS})));
NumExpts = max([OUTPUT_any.INFORMATION.experimentNum]);

DATTMP.firsttarget.PBS.FFrelbase_ExptByBin = nan(NumExpts, NumBinsToPlot); % will fill this in
DATTMP.firsttarget.MUSC.FFrelbase_ExptByBin = nan(NumExpts, NumBinsToPlot); % will fill this in
DATTMP.secondtarg.PBS.FFrelbase_ExptByBin = nan(NumExpts, NumBinsToPlot); % will fill this in
DATTMP.secondtarg.MUSC.FFrelbase_ExptByBin = nan(NumExpts, NumBinsToPlot); % will fill this in

TargList = {'firsttarget', 'secondtarg'};
DrugList = {'PBS', 'MUSC'};

for i=1:length(TargList)
    targfield = TargList{i};
    
    for ii=1:length(DrugList)
        drugfield = DrugList{ii};
        
        for j=1:NumBinsToPlot
            
            FFvals = OUTPUT_any.(targfield).FFrelBaseInBin(j).(drugfield);
            Exptnums = OUTPUT_any.INFORMATION(j).experimentNum;
           
      
            % -- STICK INTO STRUCT
            DATTMP.(targfield).(drugfield).FFrelbase_ExptByBin(Exptnums, j) = FFvals;
            
        end
    end
end

%% === PLOT (all expts)

lt_figure; hold on;

% -------------------------- FIRST TARG
targfield = 'firsttarget';
% PBS
lt_subplot(2,3,1); hold on;
title([targfield ', PBS'])
X = 1:NumBinsToPlot;
Y = DATTMP.(targfield).PBS.FFrelbase_ExptByBin;
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'b', 'LineStyle', '-'});

% MUSC
lt_subplot(2,3,2); hold on;
title([targfield ', MUSC'])
xlabel('day bin, rel consol start');
X = 1:NumBinsToPlot;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin;
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'r', 'LineStyle', '-'});

% MUSC/PBS
lt_subplot(2,3,3); hold on;
title([targfield ', MUSC/PBS'])
xlabel('day bin, rel consol start');
X = 1:NumBinsToPlot;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin./DATTMP.(targfield).PBS.FFrelbase_ExptByBin;
plot(X, Y', '-o');
ylim([-1 1]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'k', 'LineStyle', '-'});

% -------------------------- FIRST TARG
targfield = 'secondtarg';
% PBS
lt_subplot(2,3,4); hold on;
title([targfield ', PBS'])
X = 1:NumBinsToPlot;
Y = DATTMP.(targfield).PBS.FFrelbase_ExptByBin;
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'b', 'LineStyle', '-'});

% MUSC
lt_subplot(2,3,5); hold on;
title([targfield ', MUSC'])
xlabel('day bin, rel consol start');
X = 1:NumBinsToPlot;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin;
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'r', 'LineStyle', '-'});

% MUSC/PBS
lt_subplot(2,3,6); hold on;
title([targfield ', MUSC/PBS'])
xlabel('day bin, rel consol start');
X = 1:NumBinsToPlot;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin./DATTMP.(targfield).PBS.FFrelbase_ExptByBin;
plot(X, Y', '-o');
ylim([-1 1]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'k', 'LineStyle', '-'});



%% === PLOT (only expts w/full data)

lt_figure; hold on;
numbins = 3;
inds = all(~isnan(DATTMP.firsttarget.MUSC.FFrelbase_ExptByBin(:,1:numbins)),2); % good expts (rows);

% -------------------------- FIRST TARG
targfield = 'firsttarget';
% PBS
lt_subplot(2,3,1); hold on;
title([targfield ', PBS'])
X = 1:numbins;
Y = DATTMP.(targfield).PBS.FFrelbase_ExptByBin(inds, 1:numbins);
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'b', 'LineStyle', '-'});

% MUSC
lt_subplot(2,3,2); hold on;
title([targfield ', MUSC'])
xlabel('day bin, rel consol start');
X = 1:numbins;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(inds, 1:numbins);
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'r', 'LineStyle', '-'});

% MUSC/PBS
lt_subplot(2,3,3); hold on;
title([targfield ', MUSC/PBS'])
xlabel('day bin, rel consol start');
X = 1:numbins;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(inds, 1:numbins)./DATTMP.(targfield).PBS.FFrelbase_ExptByBin(inds, 1:numbins);
plot(X, Y', '-o');
ylim([-1 1]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'k', 'LineStyle', '-'});

% -------------------------- FIRST TARG
targfield = 'secondtarg';
% PBS
lt_subplot(2,3,4); hold on;
title([targfield ', PBS'])
X = 1:numbins;
Y = DATTMP.(targfield).PBS.FFrelbase_ExptByBin(inds, 1:numbins);
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'b', 'LineStyle', '-'});

% MUSC
lt_subplot(2,3,5); hold on;
title([targfield ', MUSC'])
xlabel('day bin, rel consol start');
X = 1:numbins;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(inds, 1:numbins);
plot(X, Y', '-o');
ylim([-200 250]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'r', 'LineStyle', '-'});

% MUSC/PBS
lt_subplot(2,3,6); hold on;
title([targfield ', MUSC/PBS'])
xlabel('day bin, rel consol start');
X = 1:numbins;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(inds, 1:numbins)./DATTMP.(targfield).PBS.FFrelbase_ExptByBin(inds, 1:numbins);
plot(X, Y', '-o');
ylim([-1 1]);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'k', 'LineStyle', '-'});


%% === PLOT (all expts) [summary, all on one plot]

lt_figure; hold on;
OnlyExptWithCompleteData = 0;
if OnlyExptWithCompleteData==1
numbins = 3;
inds = all(~isnan(DATTMP.firsttarget.MUSC.FFrelbase_ExptByBin(:,1:numbins)),2); % good expts (rows);
else
numbins = NumBinsToPlot;
inds = 1:NumExpts; % good expts (rows);
end

% -------------------------- FIRST TARG
targfield = 'firsttarget';
% PBS
X = 1:numbins;
Y = DATTMP.(targfield).PBS.FFrelbase_ExptByBin(inds, 1:numbins);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'b', 'LineStyle', '-'});

% MUSC
X = 1:numbins;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(inds, 1:numbins);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'r', 'LineStyle', '-'});

% -------------------------- FIRST TARG
targfield = 'secondtarg';
% PBS
X = 1:numbins;
Y = DATTMP.(targfield).PBS.FFrelbase_ExptByBin(inds, 1:numbins);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'c', 'LineStyle', '-'});

% MUSC
X = 1:numbins;
Y = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(inds, 1:numbins);
ymean = nanmean(Y,1); 
ysem = lt_sem(Y);
lt_plot(X, ymean, {'Errors', ysem, 'Color', 'm', 'LineStyle', '-'});

lt_plot_zeroline

%% === PLOT EACH DAY FOR PBS. PLOT SUMMARY FOR MUSC
OnlyKeepIfHasDataUpToLastDay = 0; % ONLY KEEPS EXPERIMENTS WITH ENOUGH DATA

lt_figure; hold on;
DaysToPlot = [0 8]; % [-2 6] means 2 days before start consol, and 6 days of consol.
BinSize = 2; % day bin.
MuscBinsToPlot = 1:4;
plotRaw = 1;
if plotRaw ==1
    YlimHz = [-250 300];
    YlimConsol = [-0.6 1.75];
else
    YlimHz = [-100 250];
    YlimConsol = [0 1.5];  
end

% --- FIRST TARG
lt_subplot(2,1,1); hold on;
title('first targ');
% Yall = nan(NumExpts, -DaysToPlot(1)+DaysToPlot(2));
Yall = [];
Yall_musc = [];
Xmusc = [];
targfield = 'firsttarget';
for exptnum=1:NumExpts
    
    FFvals = DATSTRUCT_any.(targfield)(exptnum).FFminusBase_Mean_PBS;
    consolperiod = DATSTRUCT_any.INFORMATION(exptnum).consolPeriod;
    targlearndir = DATSTRUCT_any.INFORMATION(exptnum).targ_learn_dir;
    FFvals = targlearndir*FFvals;
    
    if isempty(consolperiod)
        continue
    end
    
    inds = (consolperiod(1)+DaysToPlot(1)):(consolperiod(1)+DaysToPlot(2)-1);
    % make sure that last day does not exceed consolidation period
    if inds(end)>consolperiod(2)
        if OnlyKeepIfHasDataUpToLastDay==1
        % then not enough consol days
        continue
        end
        inds = inds(1):consolperiod(2);
    end
    
    Y = FFvals(inds);
    X = 1:length(Y);
     if plotRaw==1
   plot(X, Y, '-k');
     end
     
    % ---- combine
%     Yall(exptnum, 1:length(Y)) = Y;
    if length(Y)<DaysToPlot(2)-DaysToPlot(1)
        Y((length(Y)+1):(DaysToPlot(2)-DaysToPlot(1)))=nan;
        disp(Y);
    end
    
    Yall = [Yall; Y];


    % --- MUSC 
    FFmusc = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(exptnum, MuscBinsToPlot);
    X = (MuscBinsToPlot-1)*BinSize+(BinSize/2+0.5)-DaysToPlot(1);
      if plotRaw==1
   plot(X, FFmusc, 'o-r');
      end
      Xmusc = X;
    Yall_musc = [Yall_musc; FFmusc];
end

% mean PBS
Ymean = nanmean(Yall,1);
Ysem = lt_sem(Yall);
% lt_plot(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color','b', 'LineStyle', '-'});
shadedErrorBar(1:length(Ymean), Ymean, Ysem, {'Color','k'}, 1);
lt_plot(1:length(Ymean), Ymean, {'Color','k', 'LineStyle', '-'});


% mean MUSC
Ymean = nanmean(Yall_musc);
Ysem = lt_sem(Yall_musc);
% lt_plot(Xmusc, Ymean, {'Errors', Ysem, 'Color','r', 'LineStyle', '-'});
shadedErrorBar(Xmusc, Ymean, Ysem, {'Color','r'}, 1);
lt_plot(Xmusc, Ymean, {'Color','r', 'LineStyle', '-'});

ylim(YlimHz);
lt_plot_zeroline;
line(-[DaysToPlot(1) DaysToPlot(1)]+0.5, ylim);

% --- SECOND TARG
lt_subplot(2,1,2); hold on;
title('second targ');
% Yall = nan(NumExpts, -DaysToPlot(1)+DaysToPlot(2));
Yall = [];
Yall_musc = [];
Xmusc = [];
targfield = 'secondtarg';
for exptnum=1:NumExpts
    
    if strcmp(targfield, 'secondtarg')
    FFvals = DATSTRUCT_any.secondtarget(exptnum).FFminusBase_Mean_PBS;
    else        
    FFvals = DATSTRUCT_any.(targfield)(exptnum).FFminusBase_Mean_PBS;
    end
    consolperiod = DATSTRUCT_any.INFORMATION(exptnum).consolPeriod;
    targlearndir = DATSTRUCT_any.INFORMATION(exptnum).targ_learn_dir;
    FFvals = targlearndir*FFvals;
    
    if isempty(consolperiod)
        continue
    end
    
    inds = (consolperiod(1)+DaysToPlot(1)):(consolperiod(1)+DaysToPlot(2)-1);
    % make sure that last day does not exceed consolidation period
    if inds(end)>consolperiod(2)
        if OnlyKeepIfHasDataUpToLastDay==1
        % then not enough consol days
        continue
        end
        inds = inds(1):consolperiod(2);
    end
    
    Y = FFvals(inds);
    X = 1:length(Y);
     if plotRaw==1
    plot(X, Y, '-k');
     end
     
    % ---- combine
%     Yall(exptnum, 1:length(Y)) = Y;
    if length(Y)<DaysToPlot(2)-DaysToPlot(1)
        Y((length(Y)+1):(DaysToPlot(2)-DaysToPlot(1)))=nan;
        disp(Y);
    end
    
    Yall = [Yall; Y];


    % --- MUSC 
    FFmusc = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(exptnum, MuscBinsToPlot);
    X = (MuscBinsToPlot-1)*BinSize+(BinSize/2+0.5)-DaysToPlot(1);
      if plotRaw==1
   plot(X, FFmusc, 'o-r');
      end
      Xmusc = X;
    Yall_musc = [Yall_musc; FFmusc];
end

% mean PBS
Ymean = nanmean(Yall);
Ysem = lt_sem(Yall);
% lt_plot(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color','b', 'LineStyle', '-'});
shadedErrorBar(1:length(Ymean), Ymean, Ysem, {'Color','k'}, 1);
lt_plot(1:length(Ymean), Ymean, {'Color','k', 'LineStyle', '-'});

% mean MUSC
Ymean = nanmean(Yall_musc);
Ysem = lt_sem(Yall_musc);
% lt_plot(Xmusc, Ymean, {'Errors', Ysem, 'Color','r', 'LineStyle', '-'});
shadedErrorBar(Xmusc, Ymean, Ysem, {'Color','r'}, 1);
lt_plot(Xmusc, Ymean, {'Color','r', 'LineStyle', '-'});

ylim(YlimHz);
lt_plot_zeroline;
line(-[DaysToPlot(1) DaysToPlot(1)]+0.5, ylim);


%% === PLOT EACH DAY FOR PBS. PLOT SUMMARY FOR MUSC [NORMALIZE ALL VALUES TO MEAN PBS VAL]

lt_figure; hold on;


% --- FIRST TARG
lt_subplot(2,1,1); hold on;
title('first targ');
% Yall = nan(NumExpts, -DaysToPlot(1)+DaysToPlot(2));
Yall = [];
Yall_musc = [];
Xmusc = [];
targfield = 'firsttarget';
for exptnum=1:NumExpts
    
    FFvals = DATSTRUCT_any.(targfield)(exptnum).FFminusBase_Mean_PBS;
    consolperiod = DATSTRUCT_any.INFORMATION(exptnum).consolPeriod;
    targlearndir = DATSTRUCT_any.INFORMATION(exptnum).targ_learn_dir;
    FFvals = targlearndir*FFvals;
    
    if isempty(consolperiod)
        continue
    end
    
    inds = (consolperiod(1)+DaysToPlot(1)):(consolperiod(1)+DaysToPlot(2)-1);
    % make sure that last day does not exceed consolidation period
    if inds(end)>consolperiod(2)
        if OnlyKeepIfHasDataUpToLastDay==1
        % then not enough consol days
        continue
        end
        inds = inds(1):consolperiod(2);
    end
    
    Y = FFvals(inds);
    
    norm_Denominator = mean(Y(-DaysToPlot+1:end));  % norm factor
    
    Y = Y./norm_Denominator;
    X = 1:length(Y);
    if plotRaw==1
    plot(X, Y, '-k');
    end
    
    % ---- combine
    %     Yall(exptnum, 1:length(Y)) = Y;
    if length(Y)<DaysToPlot(2)-DaysToPlot(1)
        Y((length(Y)+1):(DaysToPlot(2)-DaysToPlot(1)))=nan;
        disp(Y);
    end
    Yall = [Yall; Y];
    
    
    % --- MUSC
    FFmusc = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(exptnum, MuscBinsToPlot);
    FFmusc = FFmusc/norm_Denominator;
    X = (MuscBinsToPlot-1)*BinSize+(BinSize/2+0.5)-DaysToPlot(1);
        if plotRaw==1
plot(X, FFmusc, 'o-r');
        end
        Xmusc = X;
    Yall_musc = [Yall_musc; FFmusc];
end

% mean PBS
Ymean = nanmean(Yall);
Ysem = lt_sem(Yall);
shadedErrorBar(1:length(Ymean), Ymean, Ysem, {'Color','k'}, 1);
lt_plot(1:length(Ymean), Ymean, {'Color','k', 'LineStyle', '-'});
% lt_plot(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color','b', 'LineStyle', '-'});

% mean MUSC
Ymean = nanmean(Yall_musc);
Ysem = lt_sem(Yall_musc);
shadedErrorBar(Xmusc, Ymean, Ysem, {'Color','r'}, 1);
lt_plot(Xmusc, Ymean, {'Color','r', 'LineStyle', '-'});
% lt_plot(Xmusc, Ymean, {'Errors', Ysem, 'Color','r', 'LineStyle', '-'});

ylim(YlimConsol);
lt_plot_zeroline;
line(-[DaysToPlot(1) DaysToPlot(1)]+0.5, ylim);


% --- SECOND TARG
lt_subplot(2,1,2); hold on;
title('second targ');
% Yall = nan(NumExpts, -DaysToPlot(1)+DaysToPlot(2));
Yall = [];
Yall_musc = [];
Xmusc = [];
targfield = 'secondtarg';
for exptnum=1:NumExpts
    
    if strcmp(targfield, 'secondtarg')
    FFvals = DATSTRUCT_any.secondtarget(exptnum).FFminusBase_Mean_PBS;
    else        
    FFvals = DATSTRUCT_any.(targfield)(exptnum).FFminusBase_Mean_PBS;
    end
    consolperiod = DATSTRUCT_any.INFORMATION(exptnum).consolPeriod;
    targlearndir = DATSTRUCT_any.INFORMATION(exptnum).targ_learn_dir;
    FFvals = targlearndir*FFvals;
    
    if isempty(consolperiod)
        continue
    end
    
    inds = (consolperiod(1)+DaysToPlot(1)):(consolperiod(1)+DaysToPlot(2)-1);
    % make sure that last day does not exceed consolidation period
    if inds(end)>consolperiod(2)
        if OnlyKeepIfHasDataUpToLastDay==1
        % then not enough consol days
        continue
        end
        inds = inds(1):consolperiod(2);
    end
    
    Y = FFvals(inds);
    
    norm_Denominator = mean(Y(-DaysToPlot+1:end));  % norm factor
    
    Y = Y./norm_Denominator;
    X = 1:length(Y);
        if plotRaw==1
plot(X, Y, 'k-');
        end
        
    % ---- combine
    %     Yall(exptnum, 1:length(Y)) = Y;
    if length(Y)<DaysToPlot(2)-DaysToPlot(1)
        Y((length(Y)+1):(DaysToPlot(2)-DaysToPlot(1)))=nan;
        disp(Y);
    end
    Yall = [Yall; Y];
    
    
    % --- MUSC
    FFmusc = DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(exptnum, MuscBinsToPlot);
    FFmusc = FFmusc/norm_Denominator;
    X = (MuscBinsToPlot-1)*BinSize+(BinSize/2+0.5)-DaysToPlot(1);
        if plotRaw==1
plot(X, FFmusc, 'o-r');
        end
        Xmusc = X;
    Yall_musc = [Yall_musc; FFmusc];
end

% mean PBS
Ymean = nanmean(Yall);
Ysem = lt_sem(Yall);
shadedErrorBar(1:length(Ymean), Ymean, Ysem, {'Color','k'}, 1);
lt_plot(1:length(Ymean), Ymean, {'Color','k', 'LineStyle', '-'});
% lt_plot(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color','b', 'LineStyle', '-'});

% mean MUSC
Ymean = nanmean(Yall_musc);
Ysem = lt_sem(Yall_musc);
shadedErrorBar(Xmusc, Ymean, Ysem, {'Color','r'}, 1);
lt_plot(Xmusc, Ymean, {'Color','r', 'LineStyle', '-'});
% lt_plot(Xmusc, Ymean, {'Errors', Ysem, 'Color','r', 'LineStyle', '-'});

ylim(YlimConsol);
lt_plot_zeroline;
line(-[DaysToPlot(1) DaysToPlot(1)]+0.5, ylim);

