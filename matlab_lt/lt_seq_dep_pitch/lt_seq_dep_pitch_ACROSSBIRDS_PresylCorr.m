function lt_seq_dep_pitch_ACROSSBIRDS_PresylAcoustic(SeqDepPitch_AcrossBirds, PARAMS, inds1)
%% ESTABLISH inds 1 - e.g. similar, diff, etc.

% inds1=SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;


%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
NumSyls=length(SeqDepPitch_AcrossBirds.AllSyllables.Learning_all);


%% Initiate plots
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


%% DO presyl tend to have higher corr ? [motif corr]
% ______________________________________ SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Motif corr');
xlabel('Pre sim ---- Pre diff'); ylabel('acoustic dist')

% ====== PLOT ACOUSTIC SIMILARITY DISTRIBUTION FOR PRESYL SIM AND DIFF
Ycorr={};
% ------------- Presyl sim
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 & inds1;

Ycorr{1}=SeqDepPitch_AcrossBirds.AllSyllables.Corr_motif_all(inds2);


% ------------ Presyl diff
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 & inds1;

Ycorr{2}=SeqDepPitch_AcrossBirds.AllSyllables.Corr_motif_all(inds2);

% -------- PLOT
lt_plot_distributionPlot(Ycorr);

distributionPlot(Ycorr,'histOpt',0, 'showMM', 6);

p=ranksum(Ycorr{1}, Ycorr{2});
lt_plot_pvalue(p, 'ranksum');



%% DO presyl tend to have higher corr ? [song corr]
% ______________________________________ SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('song corr');
xlabel('Pre sim ---- Pre diff'); ylabel('acoustic dist')

% ====== PLOT ACOUSTIC SIMILARITY DISTRIBUTION FOR PRESYL SIM AND DIFF
Ycorr={};
% ------------- Presyl sim
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 & inds1;

Ycorr{1}=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds2);


% ------------ Presyl diff
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 & inds1;

Ycorr{2}=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds2);

% -------- PLOT
lt_plot_distributionPlot(Ycorr);

distributionPlot(Ycorr,'histOpt',0, 'showMM', 6);

p=ranksum(Ycorr{1}, Ycorr{2});
lt_plot_pvalue(p, 'ranksum');


%% ============================== SCATTER PLOT OF LEARNING VS ACOUSTIC (for presyl similar and diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('corr(BK=PRESIM)');
xlabel('Song corr'); ylabel('generalization');


% ------------- Presyl sim
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 & inds1;

Ycorr=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds2);
Y_learn=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

[~,~,~,~,~,~, hplot]=lt_regress(Y_learn, Ycorr, 1, 0, 1);
set(hplot, 'Color', 'k');

% ------------ Presyl diff
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 & inds1;

Ycorr=SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all(inds2);
Y_learn=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

[~,~,~,~,~,~, hplot]=lt_regress(Y_learn, Ycorr, 1, 0, 1);
set(hplot, 'Color', 'r');


%% ============================== SCATTER PLOT OF LEARNING VS ACOUSTIC (for presyl similar and diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('corr(BK=PRESIM)');
xlabel('motif corr'); ylabel('generalization');


% ------------- Presyl sim
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 & inds1;

Ycorr=SeqDepPitch_AcrossBirds.AllSyllables.Corr_motif_all(inds2);
Y_learn=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

[~,~,~,~,~,~, hplot]=lt_regress(Y_learn, Ycorr, 1, 0, 1);
set(hplot, 'Color', 'k');

% ------------ Presyl diff
inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 & inds1;

Ycorr=SeqDepPitch_AcrossBirds.AllSyllables.Corr_motif_all(inds2);
Y_learn=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);

[~,~,~,~,~,~, hplot]=lt_regress(Y_learn, Ycorr, 1, 0, 1);
set(hplot, 'Color', 'r');






%% ======================================================== 
% % ====== PLOT ONLY RESIDUALS
% 
% % 1) ----- Perform linear regression on all data
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('linear regress all data'); ylabel('generalize'); xlabel('acoustic dist');
% Y_acoust_all=SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all(inds1);
% Y_learn_all=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds1);
% 
% [~, ~, residuals]=lt_regress(Y_learn_all, Y_acoust_all, 1, 0);
% 
% % --- one inds already used
% presyl_all=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all(inds1);
% acoustic_all=SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all(inds1);
% 
% 
% 
% % 2) ------- PLOT RESIDUALS
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Similar(BK=PRESIM) residuals');
% xlabel('acoustic dist'); ylabel('generalization');
% 
% % ------------- Presyl sim
% inds2=presyl_all==1;
% 
% Y_acoust=acoustic_all(inds2);
% Y_learn_resid=residuals(inds2);
% 
% [~,~,~,~,~,~, hplot]=lt_regress(Y_learn_resid, Y_acoust, 1);
% set(hplot, 'Color', 'k');
% 
% % ------------ Presyl diff
% inds2=presyl_all==0;
% 
% Y_acoust=acoustic_all(inds2);
% Y_learn_resid=residuals(inds2);
% 
% [~,~,~,~,~,~, hplot]=lt_regress(Y_learn_resid, Y_acoust, 1);
% set(hplot, 'Color', 'r');
% 
% 
% 
% % ===================================================================================
% % 3) ------- COMPARE RESIDUALS
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('SIMILAR residuals');
% xlabel('Pre sim ---- Pre diff'); ylabel('learn (residuals)')
% 
% Ylearn={};
% % ------------- Presyl sim
% inds2=presyl_all==1;
% 
% Ylearn{1}=residuals(inds2);
% 
% 
% % ------------ Presyl diff
% inds2=presyl_all==0;
% 
% Ylearn{2}=residuals(inds2);
% 
% 
% % -------- PLOT
% lt_plot_distributionPlot(Ylearn);
% 
% distributionPlot(Ylearn,'histOpt',0, 'showMM', 6);
% 
% p=ranksum(Ylearn{1}, Ylearn{2});
% lt_plot_pvalue(p, 'ranksum');
% 
% lt_plot_zeroline
% 
% 
% % ===================================================================================
% % 3) ------- COMPARE RESIDUALS [LIMITING TO ACOUSTIC RANGE WHERE BOTH HAVE
% % DATAPOINTS
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('residuals [limited acoustic data range]');
% xlabel('Pre sim ---- Pre diff'); ylabel('learn (residuals)')
% 
% 
% % 1) ----------------- FIGURE OUT ACOUSTIC RANGE
% Ylearn={};
% Yacoust={};
% % ------------- Presyl sim
% inds2=presyl_all==1;
% 
% Ylearn{1}=residuals(inds2);
% Yacoust{1}=acoustic_all(inds2);
% 
% % ------------ Presyl diff
% inds2=presyl_all==0;
% 
% Ylearn{2}=residuals(inds2);
% Yacoust{2}=acoustic_all(inds2);
% 
% % WHAT RANGE?
% Acoustic_range=[max([min(Yacoust{1}), min(Yacoust{2})]) min([max(Yacoust{1}), max(Yacoust{2})])];
%     disp(['acoustic range = ' num2str(Acoustic_range)]);
% 
%     
%     % 2) --------------------- GATHER DATA
% Ylearn={};
% Yacoust={};
% % ------------- Presyl sim
% inds2=presyl_all==1 & acoustic_all>Acoustic_range(1) & acoustic_all<Acoustic_range(2);
% 
% Ylearn{1}=residuals(inds2);
% Yacoust{1}=acoustic_all(inds2);
% 
% [~,~,~,~,~,~, hplot]=lt_regress(Ylearn{1}, Yacoust{1}, 1, 0, 0);
% set(hplot, 'Color', 'k');
% 
% 
% % ------------ Presyl diff
% inds2=presyl_all==0 & acoustic_all>Acoustic_range(1) & acoustic_all<Acoustic_range(2);
% 
% Ylearn{2}=residuals(inds2);
% Yacoust{2}=acoustic_all(inds2);
% 
% [~,~,~,~,~,~, hplot]=lt_regress(Ylearn{2}, Yacoust{2}, 1, 0, 0);
% set(hplot, 'Color', 'r');
% 
% 
% 
% % ------------- 4) COMPARE DISTRIBUTIONS
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('residuals [limited acoustic data range]');
% xlabel('Pre sim ---- Pre diff'); ylabel('learn (residuals)')
% 
% lt_plot_distributionPlot(Ylearn);
% 
% distributionPlot(Ylearn,'histOpt',0, 'showMM', 6);
% 
% p=ranksum(Ylearn{1}, Ylearn{2});
% lt_plot_pvalue(p, 'ranksum');
% 
% lt_plot_zeroline
% 
% %% ===== PLOT LEAERNING (not residuals)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('SIMILAR');
% xlabel('Pre sim ---- Pre diff'); ylabel('generalizatin')
% 
% Ylearn={};
% % ------------- Presyl sim
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 & inds1;
% 
% Ylearn{1}=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
% 
% % ------------ Presyl diff
% inds2=SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 & inds1;
% 
% Ylearn{2}=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all(inds2);
% 
% % -------- PLOT
% lt_plot_distributionPlot(Ylearn);
% 
% distributionPlot(Ylearn,'histOpt',0, 'showMM', 6);
% 
% p=ranksum(Ylearn{1}, Ylearn{2});
% lt_plot_pvalue(p, 'ranksum');
% 
% lt_plot_zeroline
% 

