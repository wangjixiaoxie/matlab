%% +++++++++++++++++++++++++++++++++++++++++ HISTOGRAMS OVERLAYED WITH BASELINE DRIFT



%% === all nontargets, generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(2:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'k')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'k')

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle', '--');
hbar_all=[hbar_all hbar];



% ----- LEGEND
legend(hbar_all, {'All Nontargets'})



%% === SAME TYPE generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');


% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1, 'b');
hbar_all=[hbar_all hbar];
hold on;
% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'b')

% === drift
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle','--')
hbar_all=[hbar_all hbar];



% ----- LEGEND
legend(hbar_all, {'Same Type', 'Drift [all]'})


%% === SAME TYPE SAME SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');



% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);


[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1, 'b');
hbar_all=[hbar_all hbar];
hold on;
% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'b')


% === drift
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle','--')
hbar_all=[hbar_all hbar];


% ----- LEGEND
legend(hbar_all, {'Same Type, same seq', 'Drift [all]'})


%% === SAME TYPE DIFF SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');


% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);


[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1, 'c');
hbar_all=[hbar_all hbar];
hold on;
% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'c')

% === drift
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle','--')
hbar_all=[hbar_all hbar];



% ----- LEGEND
legend(hbar_all, {'Same Type, diff seq', 'Drift [all]'})


%% === DIFF TYPE SAME SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');


% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);


[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1, 'r');
hbar_all=[hbar_all hbar];
hold on;
% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'r')

% === drift
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle','--')
hbar_all=[hbar_all hbar];




% ----- LEGEND
legend(hbar_all, {'Diff Type, same seq', 'Drift [all]'})


%% === DIFF TYPE DIFF SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');


% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);


[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1, 'm');
hbar_all=[hbar_all hbar];
hold on;
% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'm')

% === drift
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle','--')
hbar_all=[hbar_all hbar];


% ----- LEGEND
legend(hbar_all, {'Diff Type, diff seq', 'Drift [all]'})


%% === DIFF TYPE generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');


% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);


[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1, 'r');
hbar_all=[hbar_all hbar];
hold on;
% -- sample size;
lt_plot_text(max(Learning_targ_dir), 0.2, ['n=' num2str(length(Learning_targ_dir))], 'r')

% === drift
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle','--')
hbar_all=[hbar_all hbar];



% ----- LEGEND
legend(hbar_all, {'Diff Type', 'Drift [all]'})


