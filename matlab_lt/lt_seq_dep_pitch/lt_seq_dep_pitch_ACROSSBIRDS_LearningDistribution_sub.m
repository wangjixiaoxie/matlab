%% +++++++++++++++++++++++++++++++++++++++++ HISTOGRAMS OVERLAYED WITH BASELINE DRIFT



%% === all nontargets, generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'k')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
set(hbar, 'LineStyle', '--');
hbar_all=[hbar_all hbar];

% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);



% ----- LEGEND
legend(hbar_all, {'All Nontargets'})



%% === SAME TYPE generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
hbar_all=[hbar_all hbar];

% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);


% ----- LEGEND
legend(hbar_all, {'Same Type', 'Matched Drift'})


%% === SAME TYPE SAME SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% % __________________________ 1) Figure out best xbins to use
% Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
% Xcenters=Xcenters(1:2:end);
% hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, '', 1, 1, 0.2, 1);
set(hbar, 'Color', 'b')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
hbar_all=[hbar_all hbar];

% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);


% ----- LEGEND
legend(hbar_all, {'Same Type, same seq', 'Matched Drift'})


%% === SAME TYPE DIFF SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% % __________________________ 1) Figure out best xbins to use
% Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
% Xcenters=Xcenters(1:2:end);
% hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, '', 1, 1, 0.2, 1);
set(hbar, 'Color', 'c')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
hbar_all=[hbar_all hbar];


% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);

% ----- LEGEND
legend(hbar_all, {'Same Type, diff seq', 'Matched Drift'})


%% === DIFF TYPE SAME SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% % __________________________ 1) Figure out best xbins to use
% Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
% Xcenters=Xcenters(1:2:end);
% hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, '', 1, 1, 0.2, 1, 'r');
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
hbar_all=[hbar_all hbar];

% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);


% ----- LEGEND
legend(hbar_all, {'Diff Type, same seq', 'Matched Drift'})


%% === DIFF TYPE DIFF SEQ generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% % __________________________ 1) Figure out best xbins to use
% Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
% Xcenters=Xcenters(1:2:end);
% hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,Xcenters,hbar]=lt_plot_histogram(Learning_targ_dir, '', 1, 1, 0.2, 1, 'm');
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
hbar_all=[hbar_all hbar];

% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);


% ----- LEGEND
legend(hbar_all, {'Diff Type, diff seq', 'Matched Drift'})


%% === DIFF TYPE generalization
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
xlabel('Generalization [overlayed drift]');

% __________________________ 1) Figure out best xbins to use
Yall=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
[~,Xcenters,~]=lt_plot_histogram(Yall, '', 0, 0, 0);
Xcenters=Xcenters(1:2:end);
hbar_all=[];

% ___________________________ 2) Plot (COUNTS)
% ---------- all nontargs
inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
Learning_targ_dir=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds);

[~,~,hbar]=lt_plot_histogram(Learning_targ_dir, Xcenters, 1, 1, 0.2, 1);
set(hbar, 'Color', 'r')
set(hbar, 'LineStyle', '-');
hbar_all=[hbar_all hbar];

% === drift
DriftRelTarg=SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg(inds);

[~,~,hbar]=lt_plot_histogram(DriftRelTarg, Xcenters, 1, 1, 0.2, 1, 'k');
hbar_all=[hbar_all hbar];

% === text of mean, range, iqr
ymean=mean(Learning_targ_dir);
yrange=[min(Learning_targ_dir), max(Learning_targ_dir)];
yIQR=iqr(Learning_targ_dir);
lt_plot_text(1, 0.3, ['mean: ' num2str(ymean) '; range: ' num2str(yrange) '; iqr= ' num2str(yIQR)]);


% ----- LEGEND
legend(hbar_all, {'Diff Type', 'Matched Drift'})


