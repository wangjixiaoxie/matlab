function lt_seq_dep_pitch_ACROSSBIRDS_Antigen(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 7/9/16 - does difference of abs pitch at baseline predict antigeneralziation?
% i.e. anti gen is ude to use-based bias?

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);



%% === collect

isTargAll=[];
sameAll=[];
baseFFminustargAll=[];
targDirAll=[];
learnAll=[];
learnRelTargAll=[];
handlabSameAll=[];


for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        % -- targ base valu
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
            targBaseFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).meanFF_WithinTimeWindow;
            
        else
            targBaseFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).meanFF;
        end
        % --- targ base learn dir
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        
        % ===== for every nontarg collect info
        for j=1:length(SylsUnique)
            
            syl=SylsUnique{j};
            
            % --- is target?
            istarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            
            % --- same type?
            same=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            
            % --- hand lab sme?
            handlabSame=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab;
            
            % --- baseline diff from target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                baseFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
            else
                baseFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            end
            baseFFminusTarg=baseFF-targBaseFF;
            
            % --- targ shift dir
            targdir;
            
            % --- shift
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            
            % --- generalization
            learnRelTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
            
            % ============= OUTPUT
            isTargAll=[isTargAll istarg];
            sameAll=[sameAll same];
            baseFFminustargAll=[baseFFminustargAll baseFFminusTarg];
            targDirAll=[targDirAll targdir];
            learnAll=[learnAll learn];
            learnRelTargAll=[learnRelTargAll learnRelTarg];
            handlabSameAll=[handlabSameAll handlabSame];
            
        end
    end
    
end


%% ==== PLOT

% Initiate plots
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ===== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=isTargAll==0 & sameAll==1;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% ===== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=sameAll==0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');

% ======= SAMEISH TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('hand lab same');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=isTargAll==0 & handlabSameAll==1;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('hand lab diff');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=isTargAll==0 & handlabSameAll==0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');

% ================ ONLY CASES WHERE NONTARG STARTS HIGHER
% ===== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same [cases of nontarg start higher]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=isTargAll==0 & sameAll==1 & baseFFminustargAll>0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% ===== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff [cases of nontarg start higher]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=sameAll==0  & baseFFminustargAll>0;;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');

% ================ ONLY CASES WHERE NONTARG STARTS LOWER
% ===== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same [cases of nontarg start lower]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=isTargAll==0 & sameAll==1 & baseFFminustargAll<0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% ===== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff [cases of nontarg start lower]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('pitch shift, adaptive dir');

inds=sameAll==0  & baseFFminustargAll<0;;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnAll(inds).*targDirAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');



%% ==== PLOT [SAME AS ABOVE, BUT USING GENERALIZATION]

% Initiate plots
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ===== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=isTargAll==0 & sameAll==1;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% ===== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=sameAll==0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');

% ======= SAMEISH TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('hand lab same');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=isTargAll==0 & handlabSameAll==1;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('hand lab diff');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=isTargAll==0 & handlabSameAll==0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');

% ================ ONLY CASES WHERE NONTARG STARTS HIGHER
% ===== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same [cases of nontarg start higher]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=isTargAll==0 & sameAll==1 & baseFFminustargAll>0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% ===== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff [cases of nontarg start higher]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=sameAll==0  & baseFFminustargAll>0;;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');

% ================ ONLY CASES WHERE NONTARG STARTS LOWER
% ===== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same [cases of nontarg start lower]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=isTargAll==0 & sameAll==1 & baseFFminustargAll<0;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% ===== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff [cases of nontarg start lower]');
xlabel('nontarg base diff from targ (adaptive dir)');
ylabel('generalization');

inds=sameAll==0  & baseFFminustargAll<0;;
X=baseFFminustargAll(inds).*targDirAll(inds);
Y=learnRelTargAll(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');



