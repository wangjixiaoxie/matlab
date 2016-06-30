function lt_seq_dep_pitch_ACROSSBIRDS_LearningDistribution(SeqDepPitch_AcrossBirds, PARAMS)


%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
NumSyls=length(SeqDepPitch_AcrossBirds.AllSyllables.Learning_all);


%% Initiate plots
count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% Plot - histograms of z-scored learning - overlayed [COUNTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Pitch shift (z-score)');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- Targets
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'k')
hbar_all=[hbar_all hbar];

% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '--');
hbar_all=[hbar_all hbar];


% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', [0.4 0.4 0.7])
set(hbar, 'LineStyle', ':');
hbar_all=[hbar_all hbar];


% ---------- Different [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '--');
hbar_all=[hbar_all hbar];


% ---------- Different [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', [0.8 0.2 0.2])
set(hbar, 'LineStyle', ':');
hbar_all=[hbar_all hbar];


% ----- LEGEND
legend(hbar_all, {'Targets','Similar-same sequence','Similar-diff sequence','Diff-same sequence','Diff-diff sequence'})

%% PLOT - separaitng similar, diff [counts]

% ============================================= SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Pitch shift (z-score)');
title('Same-type syllables');
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- Targets
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'k')
hbar_all=[hbar_all hbar];

% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];


% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', [0.4 0.4 0.7])
set(hbar, 'LineStyle', '--');
hbar_all=[hbar_all hbar];

% ---- LEGEND
legend(hbar_all, {'Targets','Same sequence','Diff sequence'})


% ============================================ DIFFERENT
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Pitch shift (z-score)');
title('Same-type syllables');
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- Targets
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'k')
hbar_all=[hbar_all hbar];

% ---------- Different [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% ---------- Different [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', [0.8 0.2 0.2])
set(hbar, 'LineStyle', '--');
hbar_all=[hbar_all hbar];

% ---- LEGEND
legend(hbar_all, {'Targets','Same sequence','Diff sequence'})


%% PLOT - SIMILAR VS. DIFF
% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];


% ============================================= SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Pitch shift (z-score)');
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- Targets
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[Ybinned,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'k')
hbar_all=[hbar_all hbar];
% sample sizes
lt_plot_text(mean(Learning_targ_dir), 1.2*max(Ybinned), ['n=' num2str(length(Learning_targ_dir))], 'k');


% ---------- Similar
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[Ybinned,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];
% sample sizes
lt_plot_text(mean(Learning_targ_dir), 1.2*max(Ybinned), ['n=' num2str(length(Learning_targ_dir))], 'b');



% ---------- Different
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[Ybinned,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];
% sample sizes
lt_plot_text(mean(Learning_targ_dir), 1.2*max(Ybinned), ['n=' num2str(length(Learning_targ_dir))], 'r');



% ---- overlay distribution of baseline shifts
% 2.5 and 97.5 range of drift 
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line([lowerBound lowerBound], ylim, 'Color','k','LineStyle','--');
line([upperBound upperBound], ylim, 'Color','k','LineStyle','--');


% --- overlay sample sizes



% ---- LEGEND
legend(hbar_all, {'Targets','Same-type','Diff-type'})

%% Plot - histograms of z-scored learning - overlayed [PDFs]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Pitch shift (z-score)');

% ___________________________ 2) Plot (pdf)
% ---------- Targets
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'k')

% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');

% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', [0.4 0.4 0.7])
set(hbar, 'LineStyle', '--');

% ---------- Different [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');

% ---------- Different [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', [0.8 0.2 0.2])
set(hbar, 'LineStyle', '-');


%% Plot - histograms of z-scored learning relative to target [COUNTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];


% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', [0.4 0.4 0.7])
set(hbar, 'LineStyle', ':');
hbar_all=[hbar_all hbar];


% ---------- Different [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];


% ---------- Different [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', [0.8 0.2 0.2])
set(hbar, 'LineStyle', ':');
hbar_all=[hbar_all hbar];


% ----- LEGEND
legend(hbar_all, {'SameSyl-SameTransition','SameSyl-DiffTransition','DiffSyl-SameTransition','DiffSyl-DiffTransition'})




%% Plot - Genearlzation histogram (same-type, diff type)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ----- same type
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];


% ---------- Different
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];


% ----- LEGEND
legend(hbar_all, {'Same type','Diff type'})




%% === all nontargets, generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'k')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];



% ----- LEGEND
legend(hbar_all, {'All Nontargets'})


%% ============ historgram of z-score shifts
lt_figure; hold on;

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=linspace(Xcenters(1), Xcenters(end), 50);
hbar_all=[];


% ============================================= SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Pitch shift (z-score)');
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- Targets
lt_subplot(2,2,1); hold on;
title('targs');
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);
Xcenters2=Xcenters(1:2:end);
[Ybinned,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters2, 1, 0, 0.2, 1);
set(hbar, 'Color', 'k')
hbar_all=[hbar_all hbar];
% sample sizes
learnmean=mean(Learning_targ_dir);
lt_plot_text(mean(Learning_targ_dir), 1.2*max(Ybinned), ['n=' num2str(length(Learning_targ_dir)) '; mean=' num2str(learnmean)], 'k');
% ---- overlay distribution of baseline shifts
% 2.5 and 97.5 range of drift 
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line([lowerBound lowerBound], ylim, 'Color','k','LineStyle','--');
line([upperBound upperBound], ylim, 'Color','k','LineStyle','--');



% ---------- Similar
lt_subplot(2,2,2); hold on;
title('same');
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);
Xcenters2=Xcenters(1:2:end);

[Ybinned,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters2, 1, 0, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];
% sample sizes
learnmean=mean(Learning_targ_dir);
lt_plot_text(mean(Learning_targ_dir), 1.2*max(Ybinned), ['n=' num2str(length(Learning_targ_dir)) '; mean=' num2str(learnmean)], 'k');
% ---- overlay distribution of baseline shifts
% 2.5 and 97.5 range of drift 
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line([lowerBound lowerBound], ylim, 'Color','k','LineStyle','--');
line([upperBound upperBound], ylim, 'Color','k','LineStyle','--');




% ---------- Different
lt_subplot(2,2,3); hold on;
title('diff');
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;

Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);

[Ybinned,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 0, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];
% sample sizes
learnmean=mean(Learning_targ_dir);
lt_plot_text(mean(Learning_targ_dir), 1.2*max(Ybinned), ['n=' num2str(length(Learning_targ_dir)) '; mean=' num2str(learnmean)], 'k');
% ---- overlay distribution of baseline shifts
% 2.5 and 97.5 range of drift 
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line([lowerBound lowerBound], ylim, 'Color','k','LineStyle','--');
line([upperBound upperBound], ylim, 'Color','k','LineStyle','--');


%% +++++++++++++++++++++++++++++++++++++++++ HISTOGRAMS OVERLAYED WITH BASELINE DRIFT [ each own drift]

lt_seq_dep_pitch_ACROSSBIRDS_LearningDistribution_sub


%% +++++++++++++++++++++++++++++++++++++++++ HISTOGRAMS OVERLAYED WITH BASELINE DRIFT [using all nontarg drift]

lt_seq_dep_pitch_ACROSSBIRDS_LearningDistribution_sub2

%% == how many generalizers greater than or less than 1? 
greater_than_one=sum(Learning_targ_dir>=1);
less_than_negone=sum(Learning_targ_dir<=-1);

disp(['Generalizers greater than 1: ' num2str(greater_than_one)])
disp(['Generalizers less than -1: ' num2str(less_than_negone)])
disp(['Total number of syls: ' num2str(length(Learning_targ_dir))])

%% == how many generalizers greater than or less than 0.5? 
greater_than_one=sum(Learning_targ_dir>=0.5);
less_than_negone=sum(Learning_targ_dir<=-0.5);

disp(['Generalizers greater than 0.5: ' num2str(greater_than_one)])
disp(['Generalizers less than -0.5: ' num2str(less_than_negone)])
disp(['Total number of syls: ' num2str(length(Learning_targ_dir))])

%% === PLOT BARS SUMMARIZING ALL SYL TYPES [sign rank vs. 0]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('signrank');
xlabel('Generalization');

Ymeans=[];
Ysems=[];
p_versus_zero_all=[];

% -------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];

p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Diff
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Diff [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];



% ---------- Diff [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% =============== PLOT ALL
X=1:length(Ymeans);
Xlabels={'All Nontargs','SameSyl','DiffSyl','SameSyl-SameTrans','SameSyl-DiffTrans','DiffSyl-SameTrans','DiffSyl-DiffTrans'};
lt_plot_bar(X, Ymeans, {'Errors',Ysems});

% -- pvals (versus 0);
for i=1:length(p_versus_zero_all);
    p=p_versus_zero_all(i);
    
    star='';
    if p<0.05;
        star='*';
    end
    if p<0.005;
        star='**';
    end
    if p<0.0005;
        star='***';
    end
    lt_plot_text(X(i)-0.15, Ymeans(i)+0.04,star, 'r');
end
Ylim=ylim;
lt_plot_text(1, Ylim(1)+0.05,'vs. zero', 'r');

set(gca, 'XTick',X);
set(gca,'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


% ==== COMPUTE p-val comparing any two samples
% inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
% Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);
% 
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
% Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);

inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

p=ranksum(Learning_rel_targ1, Learning_rel_targ2);
disp(['p=' num2str(p)]);

%% === PLOT BARS SUMMARIZING ALL SYL TYPES [sign rank vs. 0]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('ttest');
xlabel('Generalization');

Ymeans=[];
Ysems=[];
p_versus_zero_all=[];

% -------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];

[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Diff
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Diff [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];



% ---------- Diff [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
[~, p_versus_zero]=ttest(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% =============== PLOT ALL
X=1:length(Ymeans);
Xlabels={'All Nontargs','SameSyl','DiffSyl','SameSyl-SameTrans','SameSyl-DiffTrans','DiffSyl-SameTrans','DiffSyl-DiffTrans'};
lt_plot_bar(X, Ymeans, {'Errors',Ysems});

% -- pvals (versus 0);
for i=1:length(p_versus_zero_all);
    p=p_versus_zero_all(i);
    
    star='';
    if p<0.05;
        star='*';
    end
    if p<0.005;
        star='**';
    end
    if p<0.0005;
        star='***';
    end
    lt_plot_text(X(i)-0.15, Ymeans(i)+0.04,star, 'r');
end
Ylim=ylim;
lt_plot_text(1, Ylim(1)+0.05,'vs. zero', 'r');

set(gca, 'XTick',X);
set(gca,'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


% ==== COMPUTE p-val comparing any two samples
% inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
% Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);
% 
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
% Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);

inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

p=ranksum(Learning_rel_targ1, Learning_rel_targ2);
disp(['p=' num2str(p)]);

%% === PLOT BARS SUMMARIZING ALL SYL TYPES [MODIFIED, PLOTTING SUBSET]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('generalization');

Ymeans=[];
Ysems=[];
p_versus_zero_all=[];

% % -------- all nontargs
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);
% 
% Ymeans=[Ymeans mean(Learning_rel_targ)];
% Ysems=[Ysems lt_sem(Learning_rel_targ)];
% 
% p_versus_zero=signrank(Learning_rel_targ);
% p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% % ---------- Similar
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);
% 
% Ymeans=[Ymeans mean(Learning_rel_targ)];
% Ysems=[Ysems lt_sem(Learning_rel_targ)];
% p_versus_zero=signrank(Learning_rel_targ);
% p_versus_zero_all=[p_versus_zero_all p_versus_zero];




% ---------- Similar [presim]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Similar [prediff]
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Diff
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% % ---------- Diff [presim]
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
% Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);
% 
% Ymeans=[Ymeans mean(Learning_rel_targ)];
% Ysems=[Ysems lt_sem(Learning_rel_targ)];
% p_versus_zero=signrank(Learning_rel_targ);
% p_versus_zero_all=[p_versus_zero_all p_versus_zero];
% 
% 
% 
% % ---------- Diff [prediff]
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
% Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);
% 
% Ymeans=[Ymeans mean(Learning_rel_targ)];
% Ysems=[Ysems lt_sem(Learning_rel_targ)];
% p_versus_zero=signrank(Learning_rel_targ);
% p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% =============== PLOT ALL
Colors={'b','c','r'};
Xlabels={'same-type, same-seq','same-type, diff-seq', 'diff type'};
for X=1:length(Ymeans);
lt_plot_bar(X, Ymeans(X), {'Errors',Ysems(X), 'Color', Colors{X}});
end

% -- pvals (versus 0);
for i=1:length(p_versus_zero_all);
    p=p_versus_zero_all(i);
    
    star='';
    if p<0.05;
        star='*';
    end
    if p<0.005;
        star='**';
    end
    if p<0.0005;
        star='***';
    end
    lt_plot_text(i-0.15, Ymeans(i)+0.04,star, 'r');
end
Ylim=ylim;
lt_plot_text(1, Ylim(1)+0.05,'vs. zero', 'r');

set(gca, 'XTick',1:length(Ymeans));
set(gca,'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


% ==== COMPUTE p-val comparing any two samples
inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);

inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

% inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);
% 
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
p=ranksum(Learning_rel_targ1, Learning_rel_targ2);
disp(['same(same) vs. same(diff): p=' num2str(p)]);

inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);

inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

% inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);
% 
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
p=ranksum(Learning_rel_targ1, Learning_rel_targ2);
disp(['same(diff) vs. diff: p=' num2str(p)]);





%% === PLOT BARS SUMMARIZING ALL SYL TYPES [SAME TYPE, DIFF TYPE] [GOOD]
lt_figure; hold on;
ylabel('Generalization');

Nall=[];
Ymeans=[];
Ysems=[];
p_versus_zero_all=[];


% ---------- Same
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Nall=[Nall numel(Learning_rel_targ)];
Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];


% ---------- Diff
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Nall=[Nall numel(Learning_rel_targ)];
Ymeans=[Ymeans mean(Learning_rel_targ)];
Ysems=[Ysems lt_sem(Learning_rel_targ)];
p_versus_zero=signrank(Learning_rel_targ);
p_versus_zero_all=[p_versus_zero_all p_versus_zero];




% =============== PLOT ALL
X=1:length(Ymeans);
Xlabels={'Same Type', 'Diff type'};

lt_plot_bar(1, Ymeans(1), {'Errors',Ysems(1), 'Color','b'});
lt_plot_bar(2, Ymeans(2), {'Errors',Ysems(2), 'Color','r'});

% -- plot text, with N
lt_plot_text(1, Ymeans(1)*1.1, ['N=' num2str(Nall(1)), 'mean(sem) =' num2str(Ymeans(1)) '(' num2str(Ysems(1)) ')'])
lt_plot_text(1, Ymeans(2)*1.1, ['N=' num2str(Nall(2)), 'mean(sem) =' num2str(Ymeans(2)) '(' num2str(Ysems(2)) ')'])


% -- pvals (versus 0);
for i=1:length(p_versus_zero_all);
    p=p_versus_zero_all(i);
    
    star='';
    if p<0.05;
        star='*';
    end
    if p<0.005;
        star='**';
    end
    if p<0.0005;
        star='***';
    end
    lt_plot_text(X(i)-0.15, Ymeans(i)+0.04,star, 'r');
end
Ylim=ylim;
lt_plot_text(1, Ylim(1)+0.05,'vs. zero', 'r');

set(gca, 'XTick',X);
set(gca,'XTickLabel',Xlabels);
rotateXLabels(gca, 45)

xlim([0 3])

% ==== COMPUTE p-val comparing any two samples
% inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
% Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);
% 
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
% Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
inds1=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_rel_targ1=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);

inds2=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_rel_targ2=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

p=ranksum(Learning_rel_targ1, Learning_rel_targ2);
disp(['p=' num2str(p)]);


%% ===== PLOT SEQUENCE DEPENDENCE (INCLUDING 2 BACK)
lt_figure; hold on;

Yvals={};
Ymeans=[];
Ysems=[];

% ==== same=type, both 2 back same and 1 back same
x=1;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);

% ==== same=type, 2 back same, 1-back different
x=2;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);


% ==== same=type, 2-back diff, 1 back same
x=3;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);


 % ==== same=type, 2 back and 1 back diff
x=4;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
   

 % ==== diff=type, anything back
x=5;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);



% ====== PLOT RAW
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    plot(i+0.2, Yvals{i}, 'ok');
    end
end

% ==== plot means
X=1:length(Yvals);

lt_plot_bar(X, Ymeans, {'Errors', Ysems});


% ==== stats



% ==== labels
Xlabels={'same, same','same, diff','diff, same', 'diff, diff', 'DIFF TYPE'};
set(gca, 'Xtick', X);
set(gca, 'XTickLabel', Xlabels);

rotateXLabels(gca,45)



%% ===== PLOT SEQUENCE DEPENDENCE (INCLUDING 2 BACK - IF ONEBACK DIFF, THEN DON'T CARE WHAT 2 BACK IS)
lt_figure; hold on;
lt_subplot(2,2,1); hold on;
Yvals={};
Ymeans=[];
Ysems=[];
Xlabels={'3-back','2-back','1-back', 'Diff-type'};
AcousticVals={};
CorrVals={};

% ==== same=type, both 2 back same and 1 back same
x=1;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);
Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);

acousticvals=SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all(inds);
AcousticVals{x}=acousticvals;
corrval=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds);
CorrVals=[CorrVals corrval];


% ==== same=type, 2-back diff, 1 back same
x=2;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);

acousticvals=SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all(inds);
AcousticVals{x}=acousticvals;
corrval=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds);
CorrVals=[CorrVals corrval];

 % ==== same=type, 1 back diff (don't care about 2 back)
x=3;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
   
acousticvals=SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all(inds);
AcousticVals{x}=acousticvals;
corrval=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds);
CorrVals=[CorrVals corrval];

% ==== diff=type, anything back
x=4;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);

acousticvals=SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all(inds);
AcousticVals{x}=acousticvals;
corrval=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds);
CorrVals=[CorrVals corrval];


% ====== PLOT RAW
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    plot(i+0.2, Yvals{i}, 'ok');
    end
end

% ==== plot means
X=1:length(Yvals);

lt_plot_bar(X, Ymeans, {'Errors', Ysems});
% 
% barh(X, Ymeans)

% == plot text
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    lt_plot_text(i+0.2, 1.3*max(Yvals{i}), ['n=' num2str(length(Yvals{i})) '; mean(sem)=' num2str(Ymeans(i)) '(' num2str(Ysems(i)) ')']);
    end
end


% ==== stats
% ----- difference from 0
for i=1:length(Yvals);
    if isempty(Yvals{i});
        continue;
    end
    p=signrank(Yvals{i});
    
if p<0.0005;
    lt_plot_text(i, 1.1*max(Yvals{i}), '***', 'r', 15);
elseif p<0.005;
    lt_plot_text(i, 1.1*max(Yvals{i}), '**', 'r', 15);
 elseif p<0.05;
    lt_plot_text(i, 1.1*max(Yvals{i}), '*', 'r', 15);
end
end
    

% ---- comparing pairs
for i=1:length(Yvals)
    
    for ii=1:length(Yvals);
        
        
        if i==ii;
            continue
        end
        try
            p=ranksum(Yvals{i}, Yvals{ii});
            
            if p<0.05
                disp([Xlabels{i} ' vs ' Xlabels{ii} ': p=' num2str(p)]);
            end
        catch err
        end
    end
end
      

% ---- anova (effect of sequence for same-types?)
Y=[];
Acoustic=[];
group=[];
CorrSong=[];
for i=1:3;
    Y=[Y Yvals{i}];
    Acoustic=[Acoustic AcousticVals{i}];
    group=[group i*ones(1, length(Yvals{i}))];
    CorrSong=[CorrSong CorrVals{i}];
end

% ---- linear regression (learning)
lt_subplot(2,2,2); hold on;
ylabel('gen'); xlabel('>2  back, 2back, 1back');
lt_regress(Y, group, 1, 0, 1, 1, 'k');
xlim([0 4]);

% ---- linear regression (acoustic vals)
lt_subplot(2,2,3); hold on;
ylabel('acoustic dist'); xlabel('>2  back, 2back, 1back');
lt_regress(Acoustic, group, 1, 0, 1, 1, 'k');
xlim([0 4]);

% ---- linear regression (corr vals)
lt_subplot(2,2,4); hold on;
ylabel('corr (song)'); xlabel('>2  back, 2back, 1back');
lt_regress(CorrSong, group, 1, 0, 1, 1, 'k');
xlim([0 4]);


% --- anova (learning)
keyboard
[~, ~, stats]=anovan(Y, {group});
multcompare(stats)
        
% --- anova (acoustic)
[~, ~, stats]=anovan(Acoustic, {group});
multcompare(stats)

% --- anova (corr-song)
[~, ~, stats]=anovan(CorrSong, {group});
multcompare(stats)

% --- ancova (covariate: acoustic dist)
x=Acoustic;
y=Y; % generalization
group=group;
[~, ~, ~, stats]=aoctool(x, y, group);
multcompare(stats, 'estimate', 'pmm')

   % --- ancova (covariate: corr)
x=CorrSong;
y=Y; % generalization
group=group;
[~, ~, ~, stats]=aoctool(x, y, group);
multcompare(stats, 'estimate', 'pmm')
     
           

% ==== labels
set(gca, 'Xtick', X);
set(gca, 'XTickLabel', Xlabels);
xlabel('number of syls back to first contextual diff')
rotateXLabels(gca,90)



%% ===== PLOT SEQUENCE DEPENDENCE (INCLUDING 2 BACK - IF ONEBACK DIFF, THEN DON'T CARE WHAT 2 BACK IS)
% [SEPARATING DIFF TYPES AS WELL]
lt_figure; hold on;

Yvals={};
Ymeans=[];
Ysems=[];

% ==== same=type, both 2 back same and 1 back same
x=1;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);


% ==== same=type, 2-back diff, 1 back same
x=2;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);


 % ==== same=type, 1 back diff (don't care about 2 back)
x=3;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
   

% ==== diff=type, 2 back same, 1 back same
x=4;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);

% ==== diff=type, 2 back diff, 1 back same
x=5;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);


% ==== diff=type, 1 back diff
x=6;
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);




% ====== PLOT RAW
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    plot(i+0.2, Yvals{i}, 'ok');
    end
end

% ==== plot means
X=1:length(Yvals);

lt_plot_bar(X, Ymeans, {'Errors', Ysems});
% 
% barh(X, Ymeans)

% ==== stats
% ----- difference from 0
for i=1:length(Yvals);
    if isempty(Yvals{i});
        continue;
    end
    p=signrank(Yvals{i});
    
if p<0.0005;
    lt_plot_text(i, 1.1*max(Yvals{i}), '***', 'r', 15);
elseif p<0.005;
    lt_plot_text(i, 1.1*max(Yvals{i}), '**', 'r', 15);
 elseif p<0.05;
    lt_plot_text(i, 1.1*max(Yvals{i}), '*', 'r', 15);
end
end
    

% ---- comparing pairs
for i=1:length(Yvals)
    
    for ii=1:length(Yvals);
        
        
        if i==ii;
            continue
        end
        try
            p=ranksum(Yvals{i}, Yvals{ii});
            
            if p<0.05
                disp([Xlabels{i} ' vs ' Xlabels{ii} ': p=' num2str(p)]);
            end
        catch err
        end
    end
end
        
        
           

% ==== labels
Xlabels={'3-back','2-back','1-back', '3-back', '2-back', '1-back'};
set(gca, 'Xtick', X);
set(gca, 'XTickLabel', Xlabels);
xlabel('number of syls back to first contextual diff')
rotateXLabels(gca,90)

%% ===== PLOT ALL 4 CLASSES
lt_figure; hold on;

Yvals={};
Ymeans=[];
Ysems=[];
Colors={};

% ==== same=type, same seq
x=1;
color='b';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;

% ==== same=type, pre diff
x=2;
color='c';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;


 % ==== diff type, same seq
x=3;
color='r';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
   Colors{x}=color;


% ==== diff=type, diff seq
x=4;
color='m';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;



% ====== PLOT RAW
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    plot(i+0.2, Yvals{i}, 'o', 'Color', Colors{i});
    
    lt_plot_bar(i, mean(Yvals{i}), {'Errors', lt_sem(Yvals{i}), 'Color', Colors{i}});
    end
end

% ==== plot means

% X=1:length(Yvals);
% 
% lt_plot_bar(X, Ymeans, {'Errors', Ysems});
% % 
% % barh(X, Ymeans)



% ==== labels
Xlabels={'SameType, SameTrans','SameType, DiffTrans','DiffType, SameTrans', 'DiffType, DiffTrans'};
set(gca, 'Xtick', 1:length(Xlabels));
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 90)

% ==== stats
% ----- difference from 0
for i=1:length(Yvals);
    if isempty(Yvals{i});
        continue;
    end
    [h, p]=ttest(Yvals{i});
    %      p=signrank(Yvals{i});
    disp([Xlabels{i} ' vs 0; p='  num2str(p)]);
    if p<0.0005;
        lt_plot_text(i, 1.1*max(Yvals{i}), '***', 'r', 15);
    elseif p<0.005;
        lt_plot_text(i, 1.1*max(Yvals{i}), '**', 'r', 15);
    elseif p<0.05;
        lt_plot_text(i, 1.1*max(Yvals{i}), '*', 'r', 15);
    end
end



% ---- comparing pairs
for i=1:length(Yvals)
    
    for ii=1:length(Yvals);
        
        
        if i==ii;
            continue
        end
        try
            p=ranksum(Yvals{i}, Yvals{ii});
            
            if p<0.05
                disp([Xlabels{i} ' vs ' Xlabels{ii} ': p=' num2str(p)]);
            end
        catch err
        end
    end
end
      
        

%% ===== PLOT ALL 4 CLASSES [NOT ADJACENT TO TARG]
lt_figure; hold on;

Yvals={};
Ymeans=[];
Ysems=[];
Colors={};

% ==== same=type, same seq
x=1;
color='b';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;

% ==== same=type, pre diff
x=2;
color='c';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;


 % ==== diff type, same seq
x=3;
color='r';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
   Colors{x}=color;


% ==== diff=type, diff seq
x=4;
color='m';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==0;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;



% ====== PLOT RAW
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    plot(i+0.2, Yvals{i}, 'o', 'Color', Colors{i});
    
    lt_plot_bar(i, mean(Yvals{i}), {'Errors', lt_sem(Yvals{i}), 'Color', Colors{i}});
    end
end

% ==== plot means

% X=1:length(Yvals);
% 
% lt_plot_bar(X, Ymeans, {'Errors', Ysems});
% % 
% % barh(X, Ymeans)



% ==== labels
Xlabels={'SameType, SameTrans','SameType, DiffTrans','DiffType, SameTrans', 'DiffType, DiffTrans'};
set(gca, 'Xtick', 1:length(Xlabels));
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 90)

% ==== stats
% ----- difference from 0
for i=1:length(Yvals);
    if isempty(Yvals{i});
        continue;
    end
    [h, p]=ttest(Yvals{i});
    %      p=signrank(Yvals{i});
    disp([Xlabels{i} ' vs 0; p='  num2str(p)]);
    if p<0.0005;
        lt_plot_text(i, 1.1*max(Yvals{i}), '***', 'r', 15);
    elseif p<0.005;
        lt_plot_text(i, 1.1*max(Yvals{i}), '**', 'r', 15);
    elseif p<0.05;
        lt_plot_text(i, 1.1*max(Yvals{i}), '*', 'r', 15);
    end
end



% ---- comparing pairs
for i=1:length(Yvals)
    
    for ii=1:length(Yvals);
        
        
        if i==ii;
            continue
        end
        try
            p=ranksum(Yvals{i}, Yvals{ii});
            
            if p<0.05
                disp([Xlabels{i} ' vs ' Xlabels{ii} ': p=' num2str(p)]);
            end
        catch err
        end
    end
end

title('not adjacent to targ');
        
        
           

%% ===== PLOT ALL 4 CLASSES [YES ADJACENT TO TARG]
lt_figure; hold on;

Yvals={};
Ymeans=[];
Ysems=[];
Colors={};

% ==== same=type, same seq
x=1;
color='b';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;

% ==== same=type, pre diff
x=2;
color='c';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;


 % ==== diff type, same seq
x=3;
color='r';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
   Colors{x}=color;


% ==== diff=type, diff seq
x=4;
color='m';
inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All==1;

ffvals=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

Yvals{x}=ffvals;
Ymeans(x)=mean(ffvals);
Ysems(x) = lt_sem(ffvals);
Colors{x}=color;



% ====== PLOT RAW
for i=1:length(Yvals);
    if ~isempty(Yvals{i});
        
    plot(i+0.2, Yvals{i}, 'o', 'Color', Colors{i});
    
    lt_plot_bar(i, mean(Yvals{i}), {'Errors', lt_sem(Yvals{i}), 'Color', Colors{i}});
    end
end

% ==== plot means

% X=1:length(Yvals);
% 
% lt_plot_bar(X, Ymeans, {'Errors', Ysems});
% % 
% % barh(X, Ymeans)



% ==== labels
Xlabels={'SameType, SameTrans','SameType, DiffTrans','DiffType, SameTrans', 'DiffType, DiffTrans'};
set(gca, 'Xtick', 1:length(Xlabels));
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 90)

% ==== stats
% ----- difference from 0
for i=1:length(Yvals);
    if isempty(Yvals{i});
        continue;
    end
    [h, p]=ttest(Yvals{i});
    %      p=signrank(Yvals{i});
    disp([Xlabels{i} ' vs 0; p='  num2str(p)]);
    if p<0.0005;
        lt_plot_text(i, 1.1*max(Yvals{i}), '***', 'r', 15);
    elseif p<0.005;
        lt_plot_text(i, 1.1*max(Yvals{i}), '**', 'r', 15);
    elseif p<0.05;
        lt_plot_text(i, 1.1*max(Yvals{i}), '*', 'r', 15);
    end
end



% ---- comparing pairs
for i=1:length(Yvals)
    
    for ii=1:length(Yvals);
        
        
        if i==ii;
            continue
        end
        try
            p=ranksum(Yvals{i}, Yvals{ii});
            
            if p<0.05
                disp([Xlabels{i} ' vs ' Xlabels{ii} ': p=' num2str(p)]);
            end
        catch err
        end
    end
end
        
title('adjacent to targ');

