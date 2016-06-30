function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl, epochfield_input)

%% LT 11/19/15 - Is AFP bias localized?

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var');
    norm_by_targsyl=1;
end


%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

disp('--');
disp('Removing syllables that should not be analyzed - i.e. WN overlap, since not using catch songs. REMOVED:');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(PARAMS.global.LMAN.SylsToRemove_SingleDir, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+1};
            syls_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
                    syls_unique(ind_to_remove)=[];
                    
                    % Put back into structure
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=syls_unique;
                    
                    % tell user this is done
                    disp([birdname '-' exptname ': removed ' tmp_sylremove]);
                end
            end
        end
    end
end


%% [ COLLECT DATA] - PLOT AFP BIAS VS. MP LEARNING, DISTRIBUTIONS, AND OTHER THINGS
epochfield=epochfield_input;

LearningPBS_all=[];
MPbias_all=[];
AFPbias_all=[];
SimDiff_all=[];
TargStatus_all=[];
PreSylSimilar_all=[];
Expt_count_all=[];
Yexpt_all={};
Ybird_all={};

expt_count=1;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== ONE FIGURE PER EXPT
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
            
        disp(['COLLLECTED DATA FOR : ' birdname '-' exptname]);
            
            Y_FFmean_pbs=[];
            Y_FFmean_musc=[];
            Y_FFsem_pbs=[];
            Y_FFsem_musc=[];
            Y_syls={};
            Y_similar_diff=[];
            Y_istarg=[];
            Y_AFP_bias=[];
            Y_AcousticDist=[];
            Y_Corr=[];
            Y_presimilar=[];
            Yexpt={};
            Ybird={};
            
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                
                % ===== OUTPUT DATA
                Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
%                 Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];                
%                 targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
%                 Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
                
                Yexpt=[Yexpt exptname(end-3:end)];
                Ybird=[Ybird birdname(1:4)];

                Expt_count_all=[Expt_count_all expt_count];
            end
            
            
            % ================= Flip sign if learning at targsyl is negative
            if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
                Y_FFmean_pbs=-1.*Y_FFmean_pbs;
                Y_FFmean_musc=-1.*Y_FFmean_musc;
                Y_AFP_bias=-1.*Y_AFP_bias;
            end
            
            % ========= Normalize by targsyl if desired (PBS learning
            % by taergsyl)
            if norm_by_targsyl==1;
                learning_by_targ=Y_FFmean_pbs(Y_istarg==1);
                
                Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
                Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
                Y_AFP_bias=Y_AFP_bias./learning_by_targ;
            end
            

            % ============================ COLLECT DATA TO PLOT FOR ALL
            % EXPERIMENTS
            if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
                LearningPBS_all=[LearningPBS_all Y_FFmean_pbs];
                MPbias_all=[MPbias_all Y_FFmean_musc];
                AFPbias_all=[AFPbias_all Y_AFP_bias];
                SimDiff_all=[SimDiff_all Y_similar_diff];
                TargStatus_all=[TargStatus_all Y_istarg];
                PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
                
                
                Yexpt_all=[Yexpt_all Yexpt];
                Ybird_all=[Ybird_all Ybird];

                
                expt_count=expt_count+1;
                
            end
        end
    end
end

%% ++++++++++++ PLOTS

%% ======== Regression between AFP Bias and Learning [DATA= SYLS];

count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('AFP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'Color','k'});


% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});


%% ======== Regression between AFP Bias and Learning [DATA= SYLS] [Learning>0];
learning_thresh=0;

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('AFP bias');
xlabel('learning');
title('Learning>thresh; regression');


% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'Color','k'});



% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
pvalues=[];
r2values=[];
slopevals=[];
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0 & LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'k');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================= [similar type]
% --- all nontargets
inds=TargStatus_all==0 & SimDiff_all==1 &  LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'b');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================= [DIFF]
% --- all nontargets
inds=TargStatus_all==0 & SimDiff_all==0 &  LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'r');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================ TARGS
inds=TargStatus_all==1 &  LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'y');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
lt_plot_annotation(1, {'all, similar, diff, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2)) '; ' num2str(pvalues(3)) '; ' num2str(pvalues(4))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2)) '; ' num2str(r2values(3)) '; ' num2str(r2values(4))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2)) '; ' num2str(slopevals(3)) '; ' num2str(slopevals(4))]}, 'k');



% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>learning_thresh;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});

lt_plot_zeroline



%% ======== Regression between Learning and MP bias [ALL SYLS];
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
pvalues=[];
r2values=[];
slopevals=[];
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'k');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================= [similar type]
% --- all nontargets
inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'b');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================= [DIFF]
% --- all nontargets
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'r');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================ TARGS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'y');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')

lt_plot_annotation(1, {'all, similar, diff, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2)) '; ' num2str(pvalues(3)) '; ' num2str(pvalues(4))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2)) '; ' num2str(r2values(3)) '; ' num2str(r2values(4))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2)) '; ' num2str(slopevals(3)) '; ' num2str(slopevals(4))]}, 'k');




%% ======== Regression between Learning and MP bias [ALL SYLS];
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
pvalues=[];
r2values=[];
slopevals=[];
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'k');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================= [similar type]
% --- all nontargets
inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'b');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================= [DIFF]
% --- all nontargets
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'r');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];

% ================ TARGS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'y');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')

lt_plot_annotation(1, {'all, similar, diff, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2)) '; ' num2str(pvalues(3)) '; ' num2str(pvalues(4))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2)) '; ' num2str(r2values(3)) '; ' num2str(r2values(4))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2)) '; ' num2str(slopevals(3)) '; ' num2str(slopevals(4))]}, 'k');


%% ======== Regression between Learning and MP bias [ALL SYLS] [no regression line]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')



%% ======== Regression between Learning and MP bias [ALL SYLS] [ONE REGRESSION LINE FOR ALL NONTARGS]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
pvalues=[];
r2values=[];
slopevals=[];
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'k');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ================ TARGS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'y');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')


lt_plot_annotation(1, {'all nontargs, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2))]}, 'k');



%% ======== Regression between Learning and MP bias [ALL SYLS] [ONE REGRESSION LINE FOR ALL NONTARGS] [only positive vals]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
pvalues=[];
r2values=[];
slopevals=[];
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'k');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ================ TARGS
inds=TargStatus_all==1 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'y');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')


lt_plot_annotation(1, {'all nontargs, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2))]}, 'k');



%% ======== Regression between Learning and MP bias [ALL SYLS] [REGRESSION CONSTRAINED THROUGH 0] [neg and pos values]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

tmp=X(:)\Y(:);
Ylim=ylim;
line([0 200], tmp*[0 200], 'Color','b')



% ================ TARGS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

tmp=X(:)\Y(:);
Ylim=ylim;
line([0 200], tmp*[0 200], 'Color','k')

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')


lt_plot_annotation(1, {'all nontargs, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2))]}, 'k');



%% ======== Regression between Learning and MP bias [ALL SYLS] [REGRESSION CONSTRAINED THROUGH 0] [pos values]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

tmp=X(:)\Y(:);
Ylim=ylim;
line([0 200], tmp*[0 200], 'Color','b')



% ================ TARGS
inds=TargStatus_all==1 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

tmp=X(:)\Y(:);
Ylim=ylim;
line([0 200], tmp*[0 200], 'Color','k')

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')


lt_plot_annotation(1, {'all nontargs, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2))]}, 'k');







%% ================= PLOT [LEARN(nontarg)/LEARN(targ] vs. [AFPbias(nontarg)/AFPbias(targ)]);
count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('Learn(nontarg)/Learn(targ)');
xlabel('AFPbias(nontarg)/AFPbias(targ)');

% ==== SIMILAR (presim)
inds=find(TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1);

X=[];
Y=[];
for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    afpbias=AFPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=AFPbias_all(targind);
    
    % divide nontarg by targ
    x=afpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
end

lt_plot(X, Y, {'Color','b'});



% ==== SIMILAR (prediff)
inds=find(TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0);

X=[];
Y=[];
for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    afpbias=AFPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=AFPbias_all(targind);
    
    % divide nontarg by targ
    x=afpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
end

lt_plot(X, Y, {'MarkerFaceColor', 'none', 'Color','b'});


% ==== DIFF TYPE
inds=find(TargStatus_all==0 & SimDiff_all==0);

X=[];
Y=[];
for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    afpbias=AFPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=AFPbias_all(targind);
    
    % divide nontarg by targ
    x=afpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
end

lt_plot(X, Y, {'MarkerFaceColor', 'none', 'Color','r'});



% ----------------
lt_plot_makesquare_plot45line(gca,'k')



%% ================= PLOT [LEARN(nontarg)/LEARN(targ] vs. [MPbias(nontarg)/MPbias(targ)]);
plot_text=1;
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('Learn(nontarg)/Learn(targ)');
xlabel('MPBias(nontarg)/MPBias(targ)');
title('slope<1 ==> AFP bias localized');


% ==== SIMILAR (presim)
inds=find(TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1);

X=[];
Y=[];
expt={};
bird={};
for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    mpbias=MPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=MPbias_all(targind);
    
    % divide nontarg by targ
    x=mpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
    
    expt=[expt Yexpt_all{inds(i)}];
    bird=[bird Ybird_all{inds(i)}];
   
    
end

lt_plot(X, Y, {'Color','b'});
% ----- PLOT TEXT
if plot_text==1;
    for i=1:length(X)
        
        text(X(i), Y(i), [bird{i} '-' expt{i}])
    end
end




% ==== SIMILAR (prediff)
inds=find(TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0);

X=[];
Y=[];
expt={};
bird={};

for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    mpbias=MPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=MPbias_all(targind);
    
    % divide nontarg by targ
    x=mpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
    expt=[expt Yexpt_all{inds(i)}];
    bird=[bird Ybird_all{inds(i)}];

end

lt_plot(X, Y, {'MarkerFaceColor', 'none', 'Color','b'});
% ----- PLOT TEXT
if plot_text==1;
    for i=1:length(X)
        
        text(X(i), Y(i), [bird{i} '-' expt{i}])
    end
end


% ==== DIFF TYPE
inds=find(TargStatus_all==0 & SimDiff_all==0);

X=[];
Y=[];
expt={};
bird={};

for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    mpbias=MPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=MPbias_all(targind);
    
    % divide nontarg by targ
    x=mpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
    expt=[expt Yexpt_all{inds(i)}];
    bird=[bird Ybird_all{inds(i)}];

end

lt_plot(X, Y, {'MarkerFaceColor', 'none', 'Color','r'});
% ----- PLOT TEXT
if plot_text==1;
    for i=1:length(X)
        
        text(X(i), Y(i), [bird{i} '-' expt{i}])
    end
end



% ----------------
lt_plot_makesquare_plot45line(gca,'k')



%% ================= PLOT [LEARN(nontarg)/LEARN(targ] vs. [MPbias(nontarg)/MPbias(targ)]) [SAME - no plot text]
plot_text=0;
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('Learn(nontarg)/Learn(targ)');
xlabel('MPBias(nontarg)/MPBias(targ)');
title('slope<1 ==> AFP bias localized');


% ==== SIMILAR (presim)
inds=find(TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1);

X=[];
Y=[];
expt={};
bird={};
for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    mpbias=MPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=MPbias_all(targind);
    
    % divide nontarg by targ
    x=mpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
    
    expt=[expt Yexpt_all{inds(i)}];
    bird=[bird Ybird_all{inds(i)}];
   
    
end

lt_plot(X, Y, {'Color','b'});
% ----- PLOT TEXT
if plot_text==1;
    for i=1:length(X)
        
        text(X(i), Y(i), [bird{i} '-' expt{i}])
    end
end




% ==== SIMILAR (prediff)
inds=find(TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0);

X=[];
Y=[];
expt={};
bird={};

for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    mpbias=MPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=MPbias_all(targind);
    
    % divide nontarg by targ
    x=mpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
    expt=[expt Yexpt_all{inds(i)}];
    bird=[bird Ybird_all{inds(i)}];

end

lt_plot(X, Y, {'MarkerFaceColor', 'none', 'Color','b'});
% ----- PLOT TEXT
if plot_text==1;
    for i=1:length(X)
        
        text(X(i), Y(i), [bird{i} '-' expt{i}])
    end
end


% ==== DIFF TYPE
inds=find(TargStatus_all==0 & SimDiff_all==0);

X=[];
Y=[];
expt={};
bird={};

for i=1:length(inds);
    
    learn=LearningPBS_all(inds(i));
    mpbias=MPbias_all(inds(i));
    
    % at target
    targind=TargStatus_all==1 & Expt_count_all==Expt_count_all(inds(i));
    
    learn_targ=LearningPBS_all(targind);
    afpbias_targ=MPbias_all(targind);
    
    % divide nontarg by targ
    x=mpbias/afpbias_targ;
    y=learn/learn_targ;
    
    X=[X x];
    Y=[Y y];
    expt=[expt Yexpt_all{inds(i)}];
    bird=[bird Ybird_all{inds(i)}];

end

lt_plot(X, Y, {'MarkerFaceColor', 'none', 'Color','r'});
% ----- PLOT TEXT
if plot_text==1;
    for i=1:length(X)
        
        text(X(i), Y(i), [bird{i} '-' expt{i}])
    end
end



% ----------------
lt_plot_makesquare_plot45line(gca,'k')



%% ======== EACH DAY AS DATAPOINT
%% [ COLLECT DATA] - PLOT AFP BIAS VS. MP LEARNING, DISTRIBUTIONS, AND OTHER THINGS
epochfield=epochfield_input;

LearningPBS_all=[];
MPbias_all=[];
AFPbias_all=[];
SimDiff_all=[];
TargStatus_all=[];
PreSylSimilar_all=[];
Expt_count_all=[];
DayInd_all=[];
Yexpt_all={};
Ybird_all={};

expt_count=1;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        disp(['COLLLECTED DATA FOR : ' birdname '-' exptname]);
        
        % STUFF FOR THIS EXPT
        Y_FFmean_pbs=[];
        Y_FFmean_musc=[];
        Y_FFsem_pbs=[];
        Y_FFsem_musc=[];
        Y_syls={};
        Y_similar_diff=[];
        Y_istarg=[];
        Y_AFP_bias=[];
        Y_AcousticDist=[];
        Y_Corr=[];
        Y_presimilar=[];
        Yexpt={};
        Ybird={};
        Yday=[];
        
        
        % ============ Get single-dir days with muscimol
        DaysWithMusc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        % make sure is within single dir learning
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'days_bidir_early');
            day_mustbebeforethis=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.days_bidir_early(1);

            DaysWithMusc=DaysWithMusc(DaysWithMusc<day_mustbebeforethis);
        end

        % don't take baseline
        DaysWithMusc=DaysWithMusc(DaysWithMusc>SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays(end));
        % =============
        
        
        for day = DaysWithMusc;
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            % ===== COLLECT DATA - for each syl in order, get learning (PBS and
            % MUSC)
            FF_PBS=mean(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{day}); % mean using rends across days
            FF_MUSC=mean(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day});
                        
            FFsem_PBS=lt_sem(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{day});
            FFsem_MUSC=lt_sem(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day});
            
            
            % ===== OUTPUT DATA
            Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
            Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
            Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
            Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
            Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
            Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
            Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
            Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
            Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
            %                 Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
            %                 targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            %                 Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
            
            Yexpt=[Yexpt exptname(end-3:end)];
            Ybird=[Ybird birdname(1:4)];
            Yday=[Yday day];
            Expt_count_all=[Expt_count_all expt_count];
        end
        
        end
        
        % ================= Flip sign if learning at targsyl is negative
        if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
            Y_FFmean_pbs=-1.*Y_FFmean_pbs;
            Y_FFmean_musc=-1.*Y_FFmean_musc;
            Y_AFP_bias=-1.*Y_AFP_bias;
        end
        
        % ========= Normalize by targsyl if desired (PBS learning
        % by taergsyl)
        if norm_by_targsyl==1;
            learning_by_targ=Y_FFmean_pbs(Y_istarg==1);
            
            Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
            Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
            Y_AFP_bias=Y_AFP_bias./learning_by_targ;
        end
        
        
        % ============================ COLLECT DATA TO PLOT FOR ALL
        % EXPERIMENTS
        if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
            LearningPBS_all=[LearningPBS_all Y_FFmean_pbs];
            MPbias_all=[MPbias_all Y_FFmean_musc];
            AFPbias_all=[AFPbias_all Y_AFP_bias];
            SimDiff_all=[SimDiff_all Y_similar_diff];
            TargStatus_all=[TargStatus_all Y_istarg];
            PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
            
            
            Yexpt_all=[Yexpt_all Yexpt];
            Ybird_all=[Ybird_all Ybird];
            DayInd_all=[DayInd_all  Yday];
            
            expt_count=expt_count+1;
            
        end
    end
end



%% ======== Regression between Learning and MP bias [ALL SYLS] [ONE REGRESSION LINE FOR ALL NONTARGS]
count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
ylabel('MP bias');
xlabel('learning');

% ======= TARGETS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'Color','k'});

tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT REGRESSION LINES 
pvalues=[];
r2values=[];
slopevals=[];
% =================== [all nontargs]
% --- all nontargets
inds=TargStatus_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'k');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ================ TARGS
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

[~,~,~,~,~,summarstats]=lt_regress(Y,X,1,0, 1, 0, 'y');
pvalues=[pvalues summarstats.p];
r2values=[r2values summarstats.R2];
slopevals=[slopevals summarstats.slope];


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ====== [SIMILAR PRESIM]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','b','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [SIMILAR PREDIFF]
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','b'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ====== [DIFF TYPE]
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot(X, Y, {'MarkerFaceColor','none','Color','r'});
tmp=X(:)\Y(:);
disp(['origin constrained slope: ' num2str(tmp)])


% ---------------------
lt_plot_makesquare_plot45line(gca,'k')


lt_plot_annotation(1, {'all nontargs, targ', ['p=' num2str(pvalues(1)) '; ' num2str(pvalues(2))], ...
    ['r2=' num2str(r2values(1)) '; ' num2str(r2values(2))], ...
    ['slope=' num2str(slopevals(1)) '; ' num2str(slopevals(2))]}, 'k');






