function lt_seq_dep_pitch_ACROSSBIRDS_ConsolVsBias(OUTPUTMASTER)
%% lt 5/17/17 - afp bias at 2nd targ predict consol at first targ?
% NOTE:
% bidir: uses 2nd targ (confirmed that all are saem type)
% samedir: uses 2nd targ (all expts only hjave 2 contexts, except one (rd23) for that one used the same 
%     context as in bidir (and is conservative) [confirmed all are same
%     type]
% singleTarg: mean over all Same types (different from 2nd targ)
% SUMMARY: confirmed all are same type. For samedir might want to take
% average over rd23.


%% DO NOT CHANGE 

ListOfExpts = {'OUTPUT_multidir', 'OUTPUT_samedir', 'OUTPUT_singleTarg'}; % DO NOT CHANGE

%%
ConsolidationFirstTargAll = [];
daybin_consol = 3;

MPbiasSecondTargAll = [];
MPbiasFirstTargAll = [];
daybin_mpbias = 3;

AFPbiasSecondTargAll = [];
AFPbiasFirstTargAll = [];
daybin_afpbias = 1:2;

LearningFirstTargAll = [];
LearningSecondTargAll = [];
daybin_learning = 3;

ContainsDataBinThree = [];

ExptTypeAll = []; % in order of list of expts

for k = 1:length(ListOfExpts)
    
    ExptField = ListOfExpts{k};
    OUTPUT_any = OUTPUTMASTER.(ExptField); % TEMP
    
    % =============================
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
        
        if strcmp(ExptField, 'OUTPUT_singleTarg') & strcmp(targfield, 'secondtarg');
            targfield2 = 'meanOfSameType';
        else
            targfield2 = targfield;
        end
        
        for ii=1:length(DrugList)
            drugfield = DrugList{ii};
            
            for j=1:NumBinsToPlot
                
                FFvals = OUTPUT_any.(targfield2).FFrelBaseInBin(j).(drugfield);
                Exptnums = OUTPUT_any.INFORMATION(j).experimentNum;
                
                % -- STICK INTO STRUCT
                DATTMP.(targfield).(drugfield).FFrelbase_ExptByBin(Exptnums, j) = FFvals;
                
            end
        end
    end
    
    %% ==== PLOT BARS SUMMARY (early vs. late, only keeping if has dat for both early and late)
    % averaging across days so just one val for early, one for late.
    lt_figure; hold on;
    EarlyBins = 1:2;
    LateBins = 3;
    separationdistance = 0.1;
    
    
    % ========================== FIRSTTARG
    lt_subplot(2,1,1); hold on
    targfield = 'firsttarget';
    
    % -- early
    x =1;
    
    title(targfield);
    FF_early_PBS = nanmean(DATTMP.(targfield).PBS.FFrelbase_ExptByBin(:,EarlyBins), 2);
    FF_early_MUSC = nanmean(DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(:, EarlyBins), 2);
    inds = ~isnan(nanmean(DATTMP.(targfield).PBS.FFrelbase_ExptByBin(:,LateBins), 2));
    
    X = [x-separationdistance x+separationdistance];
    Y = [FF_early_PBS(inds) FF_early_MUSC(inds)];
    plot(X, Y, '-k');
    lt_plot_bar(X, nanmean(Y,1), {'Errors', lt_sem(Y)});
    % -- put values on
    lt_plot_text(x, 1.1*max(Y(:,1)), ['mean= ' num2str(nanmean(Y,1))])
    lt_plot_text(x, max(Y(:,1)), ['sem=' num2str(lt_sem(Y))])
    % sign rank
    p = signrank(Y(:,1), Y(:,2));
    if p<0.15
        lt_plot_text(x, 1.2*max(Y(:,1)), ['p=' num2str(p, '%3.2g')], 'r');
    end
   
    
    % -- late
    x =2;
    
    FF_late_PBS = nanmean(DATTMP.(targfield).PBS.FFrelbase_ExptByBin(:,LateBins), 2);
    FF_late_MUSC = nanmean(DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(:, LateBins), 2);
    
    X = [x-separationdistance x+separationdistance];
    Y = [FF_late_PBS FF_late_MUSC];
    plot(X, Y, '-k');
    lt_plot_bar(X, nanmean(Y,1), {'Errors', lt_sem(Y)});
    % -- put values on
    lt_plot_text(x, 1.1*max(Y(:,1)), ['mean= ' num2str(nanmean(Y,1))])
    lt_plot_text(x, max(Y(:,1)), ['sem=' num2str(lt_sem(Y))])
    % sign rank
    p = signrank(Y(:,1), Y(:,2));
    if p<0.15
        lt_plot_text(x, 1.2*max(Y(:,1)), ['p=' num2str(p, '%3.2g')], 'r');
    end
    xlim([0.5 2.5])
    ylim([-200 300])
    
    % ========================== FIRSTTARG
    lt_subplot(2,1, 2); hold on
    targfield = 'secondtarg';
    
    % -- early
    x =1;
    
    title(targfield);
    FF_early_PBS = nanmean(DATTMP.(targfield).PBS.FFrelbase_ExptByBin(:,EarlyBins), 2);
    FF_early_MUSC = nanmean(DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(:, EarlyBins), 2);
    inds = ~isnan(nanmean(DATTMP.(targfield).PBS.FFrelbase_ExptByBin(:,LateBins), 2));
    
    X = [x-separationdistance x+separationdistance];
    Y = [FF_early_PBS(inds) FF_early_MUSC(inds)];
    plot(X, Y, '-k');
    lt_plot_bar(X, nanmean(Y,1), {'Errors', lt_sem(Y)});
    % -- put values on
    lt_plot_text(x, 1.1*max(Y(:,1)), ['mean= ' num2str(nanmean(Y,1))])
    lt_plot_text(x, max(Y(:,1)), ['sem=' num2str(lt_sem(Y))])
    % sign rank
    p = signrank(Y(:,1), Y(:,2));
    if p<0.15
        lt_plot_text(x, 1.2*max(Y(:,1)), ['p=' num2str(p, '%3.2g')], 'r');
    end
    
    % -- late
    x =2;
    
    FF_late_PBS = nanmean(DATTMP.(targfield).PBS.FFrelbase_ExptByBin(:,LateBins), 2);
    FF_late_MUSC = nanmean(DATTMP.(targfield).MUSC.FFrelbase_ExptByBin(:, LateBins), 2);
    
    X = [x-separationdistance x+separationdistance];
    Y = [FF_late_PBS FF_late_MUSC];
    plot(X, Y, '-k');
    lt_plot_bar(X, nanmean(Y,1), {'Errors', lt_sem(Y)});
    % -- put values on
    lt_plot_text(x, 1.1*max(Y(:,1)), ['mean= ' num2str(nanmean(Y,1))])
    lt_plot_text(x, max(Y(:,1)), ['sem=' num2str(lt_sem(Y))])
    % sign rank
    p = signrank(Y(:,1), Y(:,2));
    if p<0.15
        lt_plot_text(x, 1.2*max(Y(:,1)), ['p=' num2str(p, '%3.2g')], 'r');
    end
     xlim([0.5 2.5])
       ylim([-200 300])

    lt_subtitle(ExptField);
    
    %% ============== EXTRACT FOR CORRELATION
    consolvals = DATTMP.firsttarget.MUSC.FFrelbase_ExptByBin(:,daybin_consol)./...
        DATTMP.firsttarget.PBS.FFrelbase_ExptByBin(:,daybin_consol);
    consolvals = nanmean(consolvals, 2);
    ConsolidationFirstTargAll = [ConsolidationFirstTargAll consolvals'];
       
    mpbiasvals = DATTMP.firsttarget.MUSC.FFrelbase_ExptByBin(:, daybin_mpbias);
    mpbiasvals = nanmean(mpbiasvals, 2);
    MPbiasFirstTargAll = [MPbiasFirstTargAll mpbiasvals'];
    
    mpbiasvals = DATTMP.secondtarg.MUSC.FFrelbase_ExptByBin(:, daybin_mpbias);
    mpbiasvals = nanmean(mpbiasvals, 2);
    MPbiasSecondTargAll = [MPbiasSecondTargAll mpbiasvals'];
    
%     if strcmp(ExptField, 'OUTPUT_singleTarg')
%        fieldtmp = 'meanOfSameType'; %  
%     else
%         fieldtmp = 'secondtarg';
%     end
    afpbiasvals = DATTMP.secondtarg.PBS.FFrelbase_ExptByBin(:, daybin_afpbias) ...
        - DATTMP.secondtarg.MUSC.FFrelbase_ExptByBin(:, daybin_afpbias);
    afpbiasvals = nanmean(afpbiasvals, 2);
    AFPbiasSecondTargAll = [AFPbiasSecondTargAll afpbiasvals'];

    
    afpbiasvals = DATTMP.firsttarget.PBS.FFrelbase_ExptByBin(:, daybin_afpbias) ...
        - DATTMP.firsttarget.MUSC.FFrelbase_ExptByBin(:, daybin_afpbias);
    afpbiasvals = nanmean(afpbiasvals, 2);
    AFPbiasFirstTargAll = [AFPbiasFirstTargAll afpbiasvals'];
    
    learningtmp = DATTMP.firsttarget.PBS.FFrelbase_ExptByBin(:, daybin_learning);
    learningtmp = nanmean(learningtmp, 2);
    LearningFirstTargAll = [LearningFirstTargAll learningtmp'];
    
    learningtmp = DATTMP.secondtarg.PBS.FFrelbase_ExptByBin(:, daybin_learning);
    learningtmp = nanmean(learningtmp, 2);
    LearningSecondTargAll = [LearningSecondTargAll learningtmp'];
    
    ExptTypeAll = [ExptTypeAll k*ones(1,length(afpbiasvals))];
    
    tmp = ~isnan(DATTMP.firsttarget.PBS.FFrelbase_ExptByBin(:,3));
    ContainsDataBinThree = [ContainsDataBinThree tmp'];
    
    
    %% ====
    exptnumstmp = find(tmp==1); % those with data in bin 3
    tmpstrind = strfind(ExptField, '_');
    fieldtmp = ['DATSTRUCT_' ExptField(tmpstrind+1:end)];
    % get first day of consol (relative to start of this training);
    consolperiodstmp = [OUTPUTMASTER.(fieldtmp).INFORMATION(exptnumstmp).consolPeriod];
    ConsolFirstDays_relEpoch = consolperiodstmp(1:2:end-1);
    disp([ExptField ', consol first day (rel epoch start) = ' num2str(min(ConsolFirstDays_relEpoch)) ' to ' ...
        num2str(max(ConsolFirstDays_relEpoch))]);
    disp([ExptField ', consol first day (rel epoch start) (median) = ' num2str(median(ConsolFirstDays_relEpoch))]);
    
    
    % get first day of epoch, relative to entire trajectory
    try % since single targ doesn't have this.
    FirstDayEpoch_RelWNStart = [OUTPUTMASTER.(fieldtmp).INFORMATION(exptnumstmp).day1_FromStartWN];
    disp([ExptField ', first day of epoch rel WN start = ' num2str(min(FirstDayEpoch_RelWNStart)) ' to ' ...
        num2str(max(FirstDayEpoch_RelWNStart))]);
    disp([ExptField ', first day of epoch rel WN start (median) = ' num2str(median(FirstDayEpoch_RelWNStart))]);
    
    disp([ExptField ', N = ' num2str(length(FirstDayEpoch_RelWNStart)) ' = ' num2str(length(ConsolFirstDays_relEpoch))]);
    catch err
    end
    
    
end

%% ===== PLOT
lt_figure; hold on;

% ========== 1) raw afp bias
lt_subplot(3,2,1); hold on;
xlabel(['AFPbias, second targ (hz), bin ' num2str(daybin_afpbias)]);
ylabel(['Consolidation, bin ' num2str(daybin_consol)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll==k;
lt_plot(AFPbiasSecondTargAll(inds), ConsolidationFirstTargAll(inds), {'Color', plotcols{k}});
end
X = AFPbiasSecondTargAll;
Y = ConsolidationFirstTargAll;
lt_regress(Y, X, 0, 0, 1, 1, 'k', 1);
line(xlim, [0 0]);
line(xlim, [1 1]);
line([0 0], ylim);

legend(gca, ListOfExpts);

% ========== 2) AFP BIAS subtrast first targ AFP bias
lt_subplot(3,2,2); hold on;
xlabel(['AFPbias, second targ minus first (hz), bin ' num2str(daybin_afpbias)]);
ylabel(['Consolidation, bin ' num2str(daybin_consol)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll==k;
    X = AFPbiasSecondTargAll(inds) - AFPbiasFirstTargAll(inds);
    Y = ConsolidationFirstTargAll(inds);
    lt_plot(X, Y, {'Color', plotcols{k}});   
end
X = AFPbiasSecondTargAll - AFPbiasFirstTargAll;
Y = ConsolidationFirstTargAll;
lt_regress(Y, X, 0, 0, 1, 1, 'k', 1);
line(xlim, [0 0]);
line(xlim, [1 1]);
line([0 0], ylim);

legend(gca, ListOfExpts);

% ========== 2) (afpbias2/afopbias1)
lt_subplot(3,2,3); hold on;
xlabel(['AFPbias, secondtarg/firsttarg, bin ' num2str(daybin_afpbias)]);
ylabel(['Consolidation, bin ' num2str(daybin_consol)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll==k;
    X = AFPbiasSecondTargAll(inds)./AFPbiasFirstTargAll(inds);
    Y = ConsolidationFirstTargAll(inds);
    lt_plot(X, Y, {'Color', plotcols{k}});
end
X = AFPbiasSecondTargAll./AFPbiasFirstTargAll;
Y = ConsolidationFirstTargAll;
lt_regress(Y, X, 0, 0, 1, 1, 'k', 1);
line(xlim, [0 0]);
line(xlim, [1 1]);
line([0 0], ylim);

legend(gca, ListOfExpts);

% ========== 2) (afpbias2/afopbias1)
lt_subplot(3,2,4); hold on;
xlabel(['AFPbias, secondtarg/[abs(second)+abs(first)], bin ' num2str(daybin_afpbias)]);
ylabel(['Consolidation, bin ' num2str(daybin_consol)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll==k;
    X = AFPbiasSecondTargAll(inds)./(abs(AFPbiasSecondTargAll(inds)) + abs(AFPbiasFirstTargAll(inds)));
    Y = ConsolidationFirstTargAll(inds);
    lt_plot(X, Y, {'Color', plotcols{k}});
end
X = AFPbiasSecondTargAll./(abs(AFPbiasSecondTargAll) + abs(AFPbiasFirstTargAll));
Y = ConsolidationFirstTargAll;
lt_regress(Y, X, 0, 0, 1, 1, 'k', 1);
line(xlim, [0 0]);
line(xlim, [1 1]);
line([0 0], ylim);
line([1/2 1/2], ylim);
line([-1/2 -1/2], ylim);

legend(gca, ListOfExpts);


%% its own figure [AFP BIAS SCATTER]
lt_figure; hold on;
% ========== 2) AFP BIAS subtrast first targ AFP bias
xlabel(['AFPbias, second targ minus first (hz), bin ' num2str(daybin_afpbias)]);
ylabel(['Consolidation, bin ' num2str(daybin_consol)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll==k;
    X = AFPbiasSecondTargAll(inds) - AFPbiasFirstTargAll(inds);
    Y = ConsolidationFirstTargAll(inds);
    lt_plot(X, Y, {'Color', plotcols{k}});   
    % plot centroid
    H = errorbarxy(nanmean(X), nanmean(Y), lt_sem(X), lt_sem(Y), {'ko', plotcols{k}, plotcols{k}});
end
X = AFPbiasSecondTargAll - AFPbiasFirstTargAll;
Y = ConsolidationFirstTargAll;
lt_regress(Y, X, 0, 0, 1, 1, 'k', 1);

line(xlim, [0 0]);
line(xlim, [1 1]);
line([0 0], ylim);

legend(gca, ListOfExpts);

%% === COMPARE CONSOLIDATION BETWEEN THREE EXPERIMENTS


%% AFP BIAS FOR 3 DIFF EXPERIMENTS (BARS, PAIRED LINES)

lt_figure; hold on;

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
separation = 0.15;
title(['AFPbias(first, second), bin' num2str(daybin_afpbias)]);

for k=1:length(ListOfExpts)
    inds = ExptTypeAll==k & ContainsDataBinThree==1;
    Y = [AFPbiasFirstTargAll(inds)' AFPbiasSecondTargAll(inds)'];
    X = [k-separation k+separation];
    plot(X, Y, '-', 'Color', plotcols{k});   

    % -- plot mean
    lt_plot_bar(X, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', plotcols{k}});
    
    % -- sign rank
    p = signrank(Y(:,1), Y(:,2));
    lt_plot_text(k, max(Y(:,1)), ['p=' num2str(p)]);
end



%% ============== AFP/MP VS. GEN/SEP(learning)

% === 1) AFP CONTRIBUTION IN FIRST TARGET
Generalization = LearningSecondTargAll./LearningFirstTargAll;
AFPcontribution = AFPbiasFirstTargAll./(AFPbiasFirstTargAll + MPbiasFirstTargAll);

lt_figure; hold on;
xlabel(['Generalization, daybin ' num2str(daybin_learning)]);
ylabel(['AFPcontribution (AFP/learning) (first context), daybin ' num2str(daybin_afpbias) ' or '  num2str(daybin_mpbias)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll ==k;
    X = Generalization(inds);
    Y = AFPcontribution(inds);
    lt_plot(X, Y, {'Color', plotcols{k}})
    
    % plot centroid
    H = errorbarxy(nanmean(X), nanmean(Y), lt_sem(X), lt_sem(Y), {'ko', plotcols{k}, plotcols{k}});

end
lt_regress(AFPcontribution, Generalization, 0, 0, 1, 1, 'k', 1);

plot(xlim, [0 0]);
plot(xlim, [1 1]);
plot([0 0], ylim);
plot([1 1], ylim);

%% ============= AFP BIAS, first targ vs. second

lt_figure; hold on
ylabel(['AFP bias (2nd targ), daybin ' num2str(daybin_afpbias)]);
xlabel(['AFP bias (1st targ), daybin ' num2str(daybin_afpbias)]);

plotcols = lt_make_plot_colors(length(ListOfExpts), 0, 0);
for k=1:length(ListOfExpts)
    inds = ExptTypeAll ==k & ContainsDataBinThree==1;
    X = AFPbiasFirstTargAll(inds);
    Y = AFPbiasSecondTargAll(inds);
    lt_plot(X, Y, {'Color', plotcols{k}});
    
    % plot centroid
    H = errorbarxy(nanmean(X), nanmean(Y), lt_sem(X), lt_sem(Y), {'ko', plotcols{k}, plotcols{k}});
end
% lt_regress(AFPbiasFirstTargAll, AFPbiasSecondTargAll, 0, 0, 1, 1, 'k', 1);
xlim([-150 150]);
ylim([-150 150]);

line(xlim, [0 0]);
line([0 0], ylim);


