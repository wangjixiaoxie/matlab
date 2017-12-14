function lt_seq_dep_pitch_ACROSSBIRDS_LMANbasePlotCor(DATSTRUCT, SeqDepPitch_AcrossBirds_LMAN, ...
    PARAMS, same_type_thr, NullControl_splitday, useHandLabForSametype, ...
    onlyUseFirstExpt, recalcValues)

%% CONVERT FROM STRUCTURE TO VARIABLES (just use name of field)

fnames = fieldnames(DATSTRUCT.singlesyls);
for i=1:length(fnames)
    
   eval([fnames{i} ' = DATSTRUCT.singlesyls.(fnames{i});']);
end

fnames = fieldnames(DATSTRUCT.pairedsyls);
for i=1:length(fnames)
    
   eval([fnames{i} ' = DATSTRUCT.pairedsyls.(fnames{i});']);
end

%%
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

%% ===== MUSC drive same-type to be more correlated? [song corr]


% =========== SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

% ----- PLOT 1- SCATTER
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);

% ---- PLOT 2 - PAIRED TEST
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');

x = [1 2];
y = [CorrSong_PBS_pairs(inds)' CorrSong_MUSC_pairs(inds)'];
plot(x, y, '-ok');
plot(2, y(:,2), 'or');
% --- plot means
ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.1, ymean, {'Errors', ysem, 'LineStyle', '-', 'Color', 'y'});
% --- sign rank test
p = signrank(y(:,1), y(:,2));
lt_plot_pvalue(p, 'srank', 1);
xlim([0 3]);
ylim([-0.5 1]);
lt_plot_zeroline;


% ======================== DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

% ---------- PLOT 1 - SCATTER
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);

% -------- PLOT 2 - PAIRED TEST
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');

x = [1 2];
y = [CorrSong_PBS_pairs(inds)' CorrSong_MUSC_pairs(inds)'];
plot(x, y, '-ok');
plot(2, y(:,2), 'or');
% --- plot means
ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.1, ymean, {'Errors', ysem, 'LineStyle', '-', 'Color', 'y'});
% --- sign rank test
p = signrank(y(:,1), y(:,2));
lt_plot_pvalue(p, 'srank', 1);
xlim([0 3]);
ylim([-0.5 1]);
lt_plot_zeroline;



% =========================== COMBINED SCATTER PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME(bu), DIFF(rd)');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');

% 1) SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);


% 2) DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);


% =========================== CHANGE IN CORR RELATIVE TO STARTING CORR
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME(bu), DIFF(rd)');
xlabel('PBS - paired corr [song]');
ylabel('CHANGE IN CORR (musc - pbs)');

Xallall = [];
Yallall = [];
Grpall = [];

% 1) SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_regress(Y, X, 1, 0, 1, 1, color, 0);
lt_plot_zeroline;

Xallall = [Xallall X];
Yallall = [Yallall Y];
Grpall = [Grpall 1*ones(1,length(X))];


% 1) DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_regress(Y, X, 1, 0, 1, 1, color, 0);

Xallall = [Xallall X];
Yallall = [Yallall Y];
Grpall = [Grpall 2*ones(1,length(X))];


% ######################################### ANCOVA (change in intercept)
[h,atab,ctab,stats] = aoctool(Xallall', Yallall', Grpall', 0.05, ...
    '','','','on','parallel lines');
c = multcompare(stats)



%% ############################# difference distyributions
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('corr(MUSC)-corr(PBS)');
X=-1:0.05:1;

Yallall ={};
% same
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

Yallall{1} = Y;

% diff
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

Yallall{2} = Y;



% ================== Plot distributions side by side
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('change in corr (song)');

lt_plot_MultDist(Yallall, [1 2])
lt_plot_zeroline;

%% === figure out domain that contains bulk of data for both same and diff

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME(bu), DIFF(rd)');
xlabel('PBS - paired corr [song]');
ylabel('CHANGE IN CORR (musc - pbs)');

% 1) SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_regress(Y, X, 1, 0, 1, 1, color, 0);
lt_plot_zeroline;

% -- percentiles for this type
ptiles1 = prctile(X, [2.5 97.5]);
line([ptiles1(1) ptiles1(1)], ylim, 'Color', color);
line([ptiles1(2) ptiles1(2)], ylim, 'Color', color);


% 1) DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_regress(Y, X, 1, 0, 1, 1, color, 0);

% -- percentiles for this type
ptiles2 = prctile(X, [2.5 97.5]);
line([ptiles2(1) ptiles2(1)], ylim, 'Color', color);
line([ptiles2(2) ptiles2(2)], ylim, 'Color', color);


% ==============================
BaseCorrThreshold = [max([ptiles1(1) ptiles2(1)]) min([ptiles1(2) ptiles2(2)])];


%% ############################# difference distyributions [within overlapping domain]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('corr(MUSC)-corr(PBS)');
X=-1:0.05:1;

Yallall ={};
% same
if useHandLabForSametype==1
    inds = IsSameSyl==1 & CorrSong_PBS_pairs>BaseCorrThreshold(1) & ...
        CorrSong_PBS_pairs<BaseCorrThreshold(2);
else
    inds=AcousticDist_pairs<same_type_thr  & CorrSong_PBS_pairs>BaseCorrThreshold(1) & ...
        CorrSong_PBS_pairs<BaseCorrThreshold(2);
end
color='b';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

Yallall{1} = Y;

% diff
if useHandLabForSametype==1
    inds = IsSameSyl==0  & CorrSong_PBS_pairs>BaseCorrThreshold(1) & ...
        CorrSong_PBS_pairs<BaseCorrThreshold(2);
else
    inds=AcousticDist_pairs>=same_type_thr  & CorrSong_PBS_pairs>BaseCorrThreshold(1) & ...
        CorrSong_PBS_pairs<BaseCorrThreshold(2);
end
color='r';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

Yallall{2} = Y;



% ================== Plot distributions side by side
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('change in corr (song)');

lt_plot_MultDist(Yallall, [1 2])
lt_plot_zeroline;


%% ###################### CORR, dependence on WHETHER same motif)


% ================================ SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==1;

x = CorrSong_PBS_pairs(inds);
y = CorrSong_MUSC_pairs(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
line([-0.5 1], [-0.5 1])
% ylim([-100 100]);

aoctool(x, y, samemotif)


% ======== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff (bk = samemotif; rd = diff)');
xlabel('PBS corr');
ylabel('MUSC corr');
inds = IsSameSyl==0;

x = CorrSong_PBS_pairs(inds);
y = CorrSong_MUSC_pairs(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
line([-0.5 1], [-0.5 1])
% ylim([-100 100]);

aoctool(x, y, samemotif)


%% ======= COMPARE SAME AND DIFF ... (CORR)
x = CorrSong_PBS_pairs;
y = CorrSong_MUSC_pairs;
samesyl = IsSameSyl;

aoctool(x, y, samesyl)


%% ===================== CORR, SEPARATE BY MOTIF (DIRECTL YCOMPARE SASME AND DIFF)

ChangeInCorr = CorrSong_MUSC_pairs - CorrSong_PBS_pairs;

% ===================================== SAME MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME MOTIF (bk = same syl; rd = diff syl)');
Yall = {};

% ---- same syl
inds = IsSameMotif==1 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff syl
inds = IsSameMotif==1 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)


% ===================================== DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF MOTIF (bk = same syl; rd = diff syl)');
Yall = {};

% ---- same syl
inds = IsSameMotif==0 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff syl
inds = IsSameMotif==0 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)


%% ====================== CORR, COMPARE WITHIN OR ACROSS MOTIFS
ChangeInCorr = CorrSong_MUSC_pairs - CorrSong_PBS_pairs;

% ========== SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL (left = same motif; rt = diff motif)');
Yall = {};
ylabel('change in corr (MUSC-PBS)');

% ---- same motif
inds = IsSameMotif==1 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff motif
inds = IsSameMotif==0 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)
lt_plot_zeroline
ylim([-1 1]);
% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i-1, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(0.5, 0.8, ['(vs)p=' num2str(p)], 'b')

% ========== DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL (bk = same motif; rd = diff motif)');
Yall = {};

% ---- same motif
inds = IsSameMotif==1 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff motif
inds = IsSameMotif==0 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)
lt_plot_zeroline
ylim([-1 1]);
% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i-1, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(0.5, 0.8, ['(vs)p=' num2str(p)], 'b')

%% ######################33 CORR, DEPENDENCE ON MOTIF POSOTIION
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


ChangeInCorr = CorrSong_MUSC_pairs - CorrSong_PBS_pairs;

inds = IsSameMotif==1 & ~isnan(NumSylsInBetween);

nsyls = NumSylsInBetween(inds);
syltype = IsSameSyl(inds);
Ycorrchange = ChangeInCorr(inds);

% ========================== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
ylabel('corr (MUSC minus PBS');
xlabel('num syls in between');
x = nsyls(syltype==1);
y = Ycorrchange(syltype==1);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);
mean(y(x==0)), mean(y(x~=0));

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);

% ============================= DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');

x = nsyls(syltype==0);
y = Ycorrchange(syltype==0);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);

%% ===== MUSC drive same-type to be more correlated? [motif corr]
if (0) % REASON: focus on song corr above - can have more datapoints
    % --- 45 deg (same)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same');
    xlabel('PBS - paired corr [motif] ');
    ylabel('MUSC - paired corr [motif]');
    inds=AcousticDist_pairs<same_type_thr;
    color='b';
    
    X=CorrMotif_PBS_pairs(inds);
    Y=CorrMotif_MUSC_pairs(inds);
    lt_plot_45degScatter(X, Y, color);
    
    % --
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff');
    xlabel('PBS - paired corr [motif] ');
    ylabel('MUSC - paired corr [motif]');
    inds=AcousticDist_pairs>=same_type_thr;
    color='r';
    
    X=CorrMotif_PBS_pairs(inds);
    Y=CorrMotif_MUSC_pairs(inds);
    lt_plot_45degScatter(X, Y, color);
    
    % -----
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('corr(MUSC)-corr(PBS) [motif]');
    X=-1:0.05:1;
    % diff
    inds=AcousticDist_pairs>=same_type_thr;
    color='r';
    
    Y=CorrMotif_MUSC_pairs(inds)-CorrMotif_PBS_pairs(inds);
    lt_plot_histogram(Y, X, 1, 1, '', 1, color);
    % same
    inds=AcousticDist_pairs<same_type_thr;
    color='b';
    
    Y=CorrMotif_MUSC_pairs(inds)-CorrMotif_PBS_pairs(inds);
    lt_plot_histogram(Y, X, 1, 1, '', 1, color);
end


%% ==== effect of MUSC on corr correlate with effect of MUSC on separation?
% --- SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('corr(MUSC) - corr(PBS)[song]');
ylabel('change in abs separation');
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

% --- diff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('corr(MUSC) - corr(PBS)[song]');
ylabel('change in abs separation');
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)


disp('Effect of MUSC on pairwise corr did not predict effect of MUSC on separation [both sim and diff]');

%% ==== SAME DIRECTION VS. DIFF DIR AFP BIAS --> PREDICT CHANGE IN CORR?
% --- SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('AFP bias direction (SAME -- DIFF)');
ylabel('corr (PBS: left ... MUSC: right)');
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

Xmusc=CorrSong_MUSC_pairs(inds);
Xpbs = CorrSong_PBS_pairs(inds);
Y=DATSTRUCT.pairedsyls.similarDirOfShift_pairs(inds);

Xcell = cell(1,4);
Xcell{1} = Xpbs(Y==1); % same direction
Xcell{2} = Xmusc(Y==1); % same direction
Xcell{3} = Xpbs(Y==0); % same direction
Xcell{4} = Xmusc(Y==0); % same direction

lt_plot_MultDist(Xcell, [1 2 3 4])
lt_plot_annotation(1, 'this is opposite trend expected', 'r');

% --- diff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('AFP bias direction (SAME -- DIFF)');
ylabel('corr(MUSC) - corr(PBS)[song]');
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

Xmusc=CorrSong_MUSC_pairs(inds);
Xpbs = CorrSong_PBS_pairs(inds);
Y=DATSTRUCT.pairedsyls.similarDirOfShift_pairs(inds);

Xcell = cell(1,4);
Xcell{1} = Xpbs(Y==1); % same direction
Xcell{2} = Xmusc(Y==1); % same direction
Xcell{3} = Xpbs(Y==0); % same direction
Xcell{4} = Xmusc(Y==0); % same direction

lt_plot_MultDist(Xcell, [1 2 3 4])
% lt_plot_annotation(1, 'this is opposite trend expected', 'r');

%% === change in corr versus corr
if (0) % THIS IS REPLACED BY ABOVE (correaltion of change in corr vs. starting corr)
    figcount=1;
    subplotrows=2;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    % --
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same and diff');
    ylabel('corr(PBS) - corr(MUSC) [song]');
    xlabel('corr(MUSC)');
    inds=AcousticDist_pairs<same_type_thr;
    color='b';
    
    X=CorrSong_MUSC_pairs(inds);
    Y=CorrSong_PBS_pairs(inds)-CorrSong_MUSC_pairs(inds);
    lt_regress(Y, X, 1, 0, 1, 1, color);
    lt_plot_zeroline;
    % -
    inds=AcousticDist_pairs>=same_type_thr;
    color='r';
    
    X=CorrSong_MUSC_pairs(inds);
    Y=CorrSong_PBS_pairs(inds)-CorrSong_MUSC_pairs(inds);
    lt_regress(Y, X, 1, 0, 1, 1, color);
    lt_plot_zeroline; line([0 0], ylim, 'Color','k');
    xlim([-1 1]), ylim([-1 1]);
    
    % --
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same');
    ylabel('corr(PBS) - corr(MUSC) [song]');
    xlabel('corr(MUSC)');
    inds=AcousticDist_pairs<same_type_thr;
    color='b';
    
    X=CorrSong_MUSC_pairs(inds);
    Y=CorrSong_PBS_pairs(inds)-CorrSong_MUSC_pairs(inds);
    lt_regress(Y, X, 1, 0, 1, 1, color);
    lt_plot_zeroline; line([0 0], ylim, 'Color','k');
    xlim([-1 1]), ylim([-1 1]);
    
    % --
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff');
    ylabel('corr(PBS) - corr(MUSC) [song]');
    xlabel('corr(MUSC)');
    inds=AcousticDist_pairs>=same_type_thr;
    color='r';
    
    X=CorrSong_MUSC_pairs(inds);
    Y=CorrSong_PBS_pairs(inds)-CorrSong_MUSC_pairs(inds);
    lt_regress(Y, X, 1, 0, 1, 1, color);
    lt_plot_zeroline; line([0 0], ylim, 'Color','k');
    xlim([-1 1]), ylim([-1 1]);
end

%% ==== distributions of correlations (MUSC and PBS);

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

% -- SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';
X=-1:0.1:1;
% - PBS
Y=CorrSong_PBS_pairs(inds);
[~, ~, h]=lt_plot_histogram(Y, X, 1, 1, '', 1, color);
set(h, 'LineStyle', '--');
% - MUSC
Y=CorrSong_MUSC_pairs(inds);
[~, ~, h]=lt_plot_histogram(Y, X, 1, 1, '', 1, color);
xlim([-1 1]);


% -- DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';
X=-1:0.1:1;
% - PBS
Y=CorrSong_PBS_pairs(inds);
[~, ~, h]=lt_plot_histogram(Y, X, 1, 1, '', 1, color);
set(h, 'LineStyle', '--');
% - MUSC
Y=CorrSong_MUSC_pairs(inds);
[~, ~, h]=lt_plot_histogram(Y, X, 1, 1, '', 1, color);
xlim([-1 1]);

% --
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all nontarg');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
inds=AcousticDist_pairs>0;
color='k';
X=-1:0.1:1;
% - PBS
Y=CorrSong_PBS_pairs(inds);
[~, ~, h]=lt_plot_histogram(Y, X, 1, 1, '', 1, color);
set(h, 'LineStyle', '--');
% - MUSC
Y=CorrSong_MUSC_pairs(inds);
[~, ~, h]=lt_plot_histogram(Y, X, 1, 1, '', 1, color);
xlim([-1 1]);


%% ================== DISTRIBUTIONS, SAME PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('SAME(PBS -- MUSC) ---- DIFF (PBS -- MUSC)');
Yall = {};

% ====================== SAME
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
% - PBS
Y=CorrSong_PBS_pairs(inds);
Yall{1} = Y;
% - MUSC
Y=CorrSong_MUSC_pairs(inds);
Yall{2} = Y;


% ====================== DIFF
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
% - PBS
Y=CorrSong_PBS_pairs(inds);
Yall{3} = Y;
% - MUSC
Y=CorrSong_MUSC_pairs(inds);
Yall{4} = Y;

% ================= PLOT
lt_plot_MultDist(Yall, 0:3);