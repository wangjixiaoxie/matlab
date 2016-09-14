function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_LinModel(SeqDepPitch_AcrossBirds, PARAMS, onlyLMANexpts)
%% Hit rate does not explain pattern of genearlization?

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);



%% MAKE A NEW STRUCTURE - each datapoint is one syl (across all birds and experiments)

Learning_all=[];
Learning_TargDir_All=[];
LearningRelTarg_all=[];
Similar_all=[];

Similar_all_handlab=[];

PreSimilar_all=[];
TwoBackSimilar_all=[];
Target_all=[];
AcousticDist_all=[];
Corr_motif_all=[];
Corr_motifsubsong_all=[];
Corr_song_all=[];
BirdNum_all=[];
ExptNum_all=[];
Syl_all={};
SylSingle_all={};
OverallExptNum_all=[];
DayVals_DurWN_zscore={};
DayVals_Baseline_zscore={};
DayVals_DurWN_zscore_targdirsign={};
DayVals_Baseline_zscore_targdirsign={};

baselineDriftTargDir_All=[];
baselineDrift_RelTarg_All=[];

AdjacentToTarg_All=[];
PosRelTarg_All=[];
PosInMotif_regexp=[];
TargPosInMotif=[];

Corr_motif_all_MUSC=[];
Corr_song_all_MUSC=[];
acoustdist_MUSC_all=[];


count=0;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
                
        
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;

        % ===== ONLY LMAN?
        if onlyLMANexpts==1
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0;
                continue
            end
            
        end
            
        
        count=count+1;
        for j=1:length(syls_unique);
            syl=syls_unique{j};

            % ========== EXTRACT DATA
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            % --- Learning
            Learning_all=[Learning_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean];
            
            % learning rel target learning dir
            target_learning_sign=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            learning_targsylDir=target_learning_sign*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            Learning_TargDir_All=[Learning_TargDir_All learning_targsylDir];
            
            
            % --- Learning rel targ
            LearningRelTarg_all=[LearningRelTarg_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ];
            
            
            % ==== baseline drift
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                baselineDrift_z=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean;
                baselineDriftTargDir=target_learning_sign*baselineDrift_z;
                baselineDrift_RelTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_rel_targ;
            else
                baselineDrift_z=nan;
                baselineDriftTargDir=nan;
                baselineDrift_RelTarg=nan;
            end
          
            baselineDriftTargDir_All=[baselineDriftTargDir_All baselineDriftTargDir];
            baselineDrift_RelTarg_All=[baselineDrift_RelTarg_All baselineDrift_RelTarg];
    
            
            % --
            Similar_all = [Similar_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
            
            Similar_all_handlab= [Similar_all_handlab ...
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab];
            
            % --
            Target_all = [Target_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
            
            % ---
            AcousticDist_all = [AcousticDist_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
            
            % --- adjacent to target (pre or post)
            PosRelTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
            if PosRelTarg==-1 | PosRelTarg==1;
                AdjacentToTarg=1;
            else
                AdjacentToTarg=0;
            end
            AdjacentToTarg_All=[AdjacentToTarg_All AdjacentToTarg];
            
            PosRelTarg_All=[PosRelTarg_All PosRelTarg];
            
            % ---- pos in motif
            posinmotif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_PosInMotif_thissyl;
            PosInMotif_regexp=[PosInMotif_regexp posinmotif];
            
            targposmotif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).regexp_PosInMotif_thissyl;
            TargPosInMotif=[TargPosInMotif targposmotif];
            
            
            try
                if ~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ)
                    corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl);
                    corr_motifsub=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(targsyl);
                else
                    corr_motif= nan;
                    corr_motifsub= nan;
                end
            catch err
                disp('ERROR in corr extraction!!!');
                corr_motif= nan;
                corr_motifsub= nan;
            end
            
            Corr_motif_all=[Corr_motif_all ...
                corr_motif];
            
            Corr_motifsubsong_all=[Corr_motifsubsong_all ...
                corr_motifsub];
            
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'CORRELATIONS');
                corrsong=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            else
                corrsong=nan;
            end
            
            Corr_song_all=[Corr_song_all corrsong];
            
            % --
            BirdNum_all=[BirdNum_all i];
            
            % --
                ExptNum_all=[ExptNum_all ii];
                
                %
                Syl_all=[Syl_all syl];
                
                % single syl
                SylSingle_all=[SylSingle_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
               
                
                % preceding syllable
                presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                PreSimilar_all=[PreSimilar_all presimilar];
                
                
                % two back
                TwoBackSim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ;
                TwoBackSimilar_all=[TwoBackSimilar_all TwoBackSim];

                
                % 
                OverallExptNum_all=[OverallExptNum_all count];
            
                % day by day pitch (zscore) vals - during WN
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                WNoffInd=min([length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore), ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd]); % min([length(data) WNoffInd]);
                dayvals1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore(WNonInd:WNoffInd);
                
                DayVals_DurWN_zscore=[DayVals_DurWN_zscore dayvals1];
                
                % day by day pitch (zscore) vals - during baseline
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                dayvals2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore(1:WNonInd-1);
                
                DayVals_Baseline_zscore=[DayVals_Baseline_zscore dayvals2];
                
                % day by day pitch (zscore) vals - during WN [learning dir
                % sign]
                DayVals_DurWN_zscore_targdirsign=[DayVals_DurWN_zscore_targdirsign target_learning_sign*dayvals1];
                
                % day by day pitch (zscore) vals - during baseline
                % [learning dir sign]
                DayVals_Baseline_zscore_targdirsign=[DayVals_Baseline_zscore_targdirsign target_learning_sign*dayvals2];
                
                % ===== LMAN STATS (if not LMAN expt, then give a nan)
                % - corr
                corr_motif_subsong_musc=nan;
                try
                corr_motif_subsong_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(targsyl);
                catch err
                end
                Corr_motif_all_MUSC=[Corr_motif_all_MUSC corr_motif_subsong_musc];
                
                corr_song_musc=nan;
                try
                corr_song_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
                catch err
                end
                Corr_song_all_MUSC=[Corr_song_all_MUSC corr_song_musc];
                
                % -- acoustic
                acoustdist_MUSC=nan;
                try
                    acoustdist_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).eucldist_from_targ_zscore;
                catch err
                end
                acoustdist_MUSC_all=[acoustdist_MUSC_all acoustdist_MUSC];
                
                
                
                
                
        end
    end
end

% ----- OUTPUT
AllSylsStruct.Learning_all=Learning_all;
AllSylsStruct.Learning_TargDir_All=Learning_TargDir_All;
AllSylsStruct.LearningRelTarg_all=LearningRelTarg_all;
AllSylsStruct.Similar_all=Similar_all;
AllSylsStruct.Similar_all_handlab=Similar_all_handlab;
AllSylsStruct.PreSimilar_all=PreSimilar_all;
AllSylsStruct.TwoBackSimilar_all=TwoBackSimilar_all;
AllSylsStruct.Target_all=Target_all;
AllSylsStruct.AcousticDist_all=AcousticDist_all;
AllSylsStruct.BirdNum_all=BirdNum_all;
AllSylsStruct.ExptNum_all=ExptNum_all;
AllSylsStruct.Syl_all=Syl_all;
AllSylsStruct.SylSingle_all=SylSingle_all;
AllSylsStruct.OverallExptNum_all=OverallExptNum_all;
AllSylsStruct.Corr_motif_all=Corr_motif_all;
AllSylsStruct.Corr_motifsubsong_all=Corr_motifsubsong_all;
AllSylsStruct.Corr_song_all=Corr_song_all;

AllSylsStruct.AdjacentToTarg_All=AdjacentToTarg_All;
AllSylsStruct.PosRelTarg_All=PosRelTarg_All;
AllSylsStruct.PosInMotif_regexp=PosInMotif_regexp;
AllSylsStruct.TargPosInMotif=TargPosInMotif;

AllSylsStruct.baselineDriftTargDir_Z=baselineDriftTargDir_All;
AllSylsStruct.baselineDriftRelTarg=baselineDrift_RelTarg_All;
AllSylsStruct.DayVals_DurWN_zscore=DayVals_DurWN_zscore;
AllSylsStruct.DayVals_Baseline_zscore=DayVals_Baseline_zscore;
AllSylsStruct.DayVals_DurWN_zscore_targdirsign=DayVals_DurWN_zscore_targdirsign;
AllSylsStruct.DayVals_Baseline_zscore_targdirsign=DayVals_Baseline_zscore_targdirsign;

AllSylsStruct.LMAN.Corr_motif_all_MUSC=Corr_motif_all_MUSC;
AllSylsStruct.LMAN.Corr_song_all_MUSC=Corr_song_all_MUSC;
AllSylsStruct.LMAN.AcousticDist_all=acoustdist_MUSC_all;


%% ++++++++++++++++
%% ====== PRESIM (effect independent of corr and acoustic and motif?)

% === plot bars
lt_figure; hold on;
title('same(presim,predif) - diff(presim, prediff)');
ylabel('generalization');

% -- similar (presim)
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & AllSylsStruct.PreSimilar_all==1;
color='b';
x=1;

y=AllSylsStruct.LearningRelTarg_all(inds);
plot(x+0.1, y, '.', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y), 'Color', color});
p=signrank(y);
if p<0.1
    lt_plot_text(x, mean(y)+0.1, ['p=' num2str(p, '%3.2g')]);
end

% -- similar (predif)
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & AllSylsStruct.PreSimilar_all==0;
color='c';
x=2;

y=AllSylsStruct.LearningRelTarg_all(inds);
plot(x+0.1, y, '.', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y), 'Color', color});
p=signrank(y);
if p<0.1
    lt_plot_text(x, mean(y)+0.1, ['p=' num2str(p, '%3.2g')]);
end
% -- diff (presim)
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & AllSylsStruct.PreSimilar_all==1;
color='r';
x=3;

y=AllSylsStruct.LearningRelTarg_all(inds);
plot(x+0.1, y, '.', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y), 'Color', color});
p=signrank(y);
if p<0.1
    lt_plot_text(x, mean(y)+0.1, ['p=' num2str(p, '%3.2g')]);
end
% -- diff (prediff)
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & AllSylsStruct.PreSimilar_all==0;
color='m';
x=4;

y=AllSylsStruct.LearningRelTarg_all(inds);
plot(x+0.1, y, '.', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y), 'Color', color});
p=signrank(y);
if p<0.1
    lt_plot_text(x, mean(y)+0.1, ['p=' num2str(p, '%3.2g')]);
end



% --- account for acoustic
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.PreSimilar_all(inds);
Zname='presim?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');



% --- account for corr
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.Corr_song_all(inds);
Xname='corr(song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.PreSimilar_all(inds);
Zname='presim?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');


% --- account for motif (same or diff)
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=~isnan(AllSylsStruct.PosRelTarg_All(inds)); % 1 = on same motif
Xname='same motif?';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.PreSimilar_all(inds);
Zname='presim?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');


%% ACOUSTIC DIST
% ======== ANCOVA (one way) (gener vs. acoustic)

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');



% === linear regression [learning vs. acoustic dist]
lt_figure; hold on;
% ---
lt_subplot(2,2,1); hold on
ylabel('gener');
xlabel('acoustic dist');
title('same-type');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

lt_regress(Y, X, 1, 0, 1, 1, 'k');

% ---
lt_subplot(2,2,2); hold on
ylabel('gener');
xlabel('acoustic dist');
title('diff-type');

inds=AllSylsStruct.Similar_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

lt_regress(Y, X, 1, 0, 1, 1, 'k');

% ==
lt_subplot(2,2,3); hold on
ylabel('gener');
xlabel('acoustic dist');
title('all nontarg');

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

lt_regress(Y, X, 1, 0, 1, 1, 'k');


disp('From ANCOVA looks like slope for same-type is greater than for diff type');
disp('Hoever, from linear regression neither slope is significant (why? - from ancova it is different even from 0)');


%% +++ ACOUSTIC DIST [SAME AS ABOVE, BUT USING HAND LABELED]

% === linear regression [learning vs. acoustic dist]
lt_figure; hold on;

% --- SAME
lt_subplot(2,2,1); hold on
ylabel('gener');
xlabel('acoustic dist');
title('same-type [hand lab]');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all_handlab==1;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

lt_regress(Y, X, 1, 0, 1, 1, 'k');


% --- DIFF
lt_subplot(2,2,2); hold on
ylabel('gener');
xlabel('acoustic dist');
title('diff-type [hand lab]');

inds=AllSylsStruct.Similar_all_handlab==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

lt_regress(Y, X, 1, 0, 1, 1, 'k');


% == ALL
lt_subplot(2,2,3); hold on
ylabel('gener');
xlabel('acoustic dist');
title('all nontarg [hand lab]');

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

lt_regress(Y, X, 1, 0, 1, 1, 'k');





%% +++++++++++++++++++++++++ correlations
% ===================== ANCOVA (one way) (gener vs. corr (song))

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.Corr_song_all(inds);
Xname='corr (song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);

disp('slope for diff type is sign. positive? for same-type is no effect');


% ===== linear regression
lt_figure; hold on;
% --- same type
lt_subplot(2,2,1); hold on;
title('same type');
xlabel('corr(song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');
lt_plot_zeroline;
lt_plot_zeroline_vert;


% -- diff type
lt_subplot(2,2,2); hold on;
title('diff type');
xlabel('corr (song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');
lt_plot_zeroline;
lt_plot_zeroline_vert;



% --- all combined
lt_subplot(2,2,3); hold on;
title('all nontarg');
xlabel('corr (song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ===== linear regression [same motif]
lt_figure; hold on;
% --- same type
lt_subplot(2,2,1); hold on;
title('same type');
xlabel('corr(motif sub song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');

% -- diff type
lt_subplot(2,2,2); hold on;
title('diff type');
xlabel('corr (motif sub song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');

% -- all nontarg
lt_subplot(2,2,3); hold on;
title('all nontarg');
xlabel('corr (motif sub song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');
line(xlim, [0 0]);
line([0 0], ylim);



% ===== linear regression [diff motif]
lt_figure; hold on;
% --- same type
lt_subplot(2,2,1); hold on;
title('same');
xlabel('corr(song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');
% --- diff type
lt_subplot(2,2,2); hold on;
title('diff');
xlabel('corr(song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');
% --- all type
lt_subplot(2,2,3); hold on;
title('all nontarg');
xlabel('corr(song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');




% ============ ANCOVA (one way) (gener vs. corr (same motif motif))

inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.Corr_motif_all);

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Xname='corr (motif)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');


% ============ ANCOVA (one way) (same motif, using song corr)

inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.Corr_motif_all);

X=AllSylsStruct.Corr_song_all(inds);
Xname='corr (song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');

% ======= ANCOVA (one way) (gener vs. corr (diff motif))

inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.Corr_song_all(inds);
Xname='corr (song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);


% ============= is the slope for diff type types (learning vs. corr) different depending on if on same or diff motif?
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;


X=AllSylsStruct.Corr_song_all(inds);
Xname='corr (song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=~isnan(AllSylsStruct.PosRelTarg_All(inds));
Zname='same-motif?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');


% ============= is the slope for same types (learning vs. corr) different depending on if on same or diff motif?
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all;


X=AllSylsStruct.Corr_song_all(inds);
Xname='corr (song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=~isnan(AllSylsStruct.PosRelTarg_All(inds));
Zname='same-motif?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');


% ============= is the slope for same types [diff presim] (learning vs. corr) different depending on if on same or diff motif?
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all & AllSylsStruct.PreSimilar_all==0;

X=AllSylsStruct.Corr_song_all(inds);
Xname='corr (song)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=~isnan(AllSylsStruct.PosRelTarg_All(inds));
Zname='same-motif?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
try
[c, m, h, nms]=multcompare(stats, 'estimate', 'slope');
catch 
    err
end

disp('===== CONCLUSIONS (trial by trial corr)');
disp('positive correlation between learning and corr for same motif (dff type only), but not for diff motif (even if use song corr for both analyses) (even though same-types do have greater mean corr for diff motif)');
disp('however, there is no sign diff between slopes (for same-types or diff types) of learn vs. corr depending on motif');
disp('note, howeber, that slopes of same-type are not significantly different (close though) than diff type (for same-motif)');


%% ==== PITCH CORREALTION [AS ABOVE, BUT USING HAND LABELED]

% ===== linear regression
lt_figure; hold on;
% --- same type
lt_subplot(2,2,1); hold on;
title('same type [HAND LAB]');
xlabel('corr(song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all_handlab==1;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% -- diff type
lt_subplot(2,2,2); hold on;
title('diff type [HAND LAB]');
xlabel('corr (song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all_handlab==0;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% --- all combined
lt_subplot(2,2,3); hold on;
title('all nontarg [HAND LAB]');
xlabel('corr (song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');

lt_plot_zeroline;
lt_plot_zeroline_vert;



% ========================== SAME MOTIF
lt_figure; hold on;
% --- same type
lt_subplot(2,2,1); hold on;
title('same type [HAND LAB]');
xlabel('corr(motif sub song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all_handlab==1;

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% -- diff type
lt_subplot(2,2,2); hold on;
title('diff type [HAND LAB]');
xlabel('corr (motif sub song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all_handlab==0;

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% -- all nontarg
lt_subplot(2,2,3); hold on;
title('all nontarg [HAND LAB]');
xlabel('corr (motif sub song)');
ylabel('gen');

inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.Corr_motifsubsong_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k');
lt_plot_zeroline;
lt_plot_zeroline_vert;



%% == indepnednet contriubtions by corr and acoustic to generalization?
lt_figure; hold on;
% -- all nontargs
lt_subplot(2,2,1); hold on;
title('all nontarg');
xlabel('corr (song)');
ylabel('acoust dist');
zlabel('generalization');
inds=AllSylsStruct.Target_all==0;

X1=AllSylsStruct.Corr_song_all(inds)'; % corr
X2=AllSylsStruct.AcousticDist_all(inds)'; % acoustic
Y=AllSylsStruct.LearningRelTarg_all(inds)'; % gener

% XX=[ones(size(X1)) X1 X2 X1.*X2];
XX=[ones(size(X1)) X1 X2];
[b, bint, ~, ~, stats]= regress(Y, XX);

scatter3(X1, X2, Y, 'filled');
x1fit = linspace(min(X1), max(X1), 50);
x2fit = linspace(min(X1), max(X2), 50);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
grid on;

lt_plot_annotation(1, ['R2=' num2str(stats(1)) '; b(corr) CI = '  num2str(bint(2,:)) '; b(acoust) CI = '  num2str(bint(3,:))], 'k')

% =========================== correlation between corr and acoustic?
lt_figure; hold on;
% -- all nontargs
lt_subplot(2,2,2); hold on;
title('all nontarg');
xlabel('corr (song)');
ylabel('acoust dist');
inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.Corr_song_all(inds)'; % corr
Y=AllSylsStruct.AcousticDist_all(inds)'; % acoustic
lt_regress(Y, X, 1, 0, 1, 1, 'k');




%% ++++++++++++++++++++ MOTIF EFFECT [look at conclusion below]

% ========= plot motif position
lt_figure; hold on;
% ==== same type
% 1) summary
lt_subplot(3,1,1); hold on;
title('same type'); xlabel('pos rel targ'); ylabel('gener');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & ~isnan(AllSylsStruct.PosRelTarg_All);
color='b';

X=AllSylsStruct.PosRelTarg_All(inds); % pos rel targ
Y=AllSylsStruct.LearningRelTarg_all(inds);
plot(X, Y, 'o', 'Color',color);
Xlim=xlim;
lt_plot_zeroline;
set(gca, 'XTick', [min(X):max(X)])
% -anova
p=anovan(Y, {X});
lt_plot_pvalue(p, 'anova', 1);
% - plot mean
distanceslist=unique(X);
distanceslist=distanceslist(~isnan(distanceslist));
for i=1:length(distanceslist);
    dist=distanceslist(i);
    
    ymean=mean(Y(X==dist));
    ysem=lt_sem(Y(X==dist));
    lt_plot_bar(dist, ymean, {'Errors', ysem, 'Color',color});
end
% 2) linear regression (hoff style) (preceding)
lt_subplot(3,2,3); hold on;
xlabel('position'); ylabel('generalization');
title('preceeding targ');
X2=X(X<0);
Y2=Y(X<0);
lt_regress(Y2, X2, 1, 0, 1, 1, 'b');
xlim([Xlim(1)-1 0]); lt_plot_zeroline;
% 3) linear regression (hoff style) (following)
lt_subplot(3,2,4); hold on;
title('following targ');
X2=X(X>0);
Y2=Y(X>0);
lt_regress(Y2, X2, 1, 0, 1, 1, 'b');
xlim([0 Xlim(2)+1]); lt_plot_zeroline;
% 4) learning vs. absolute pos from targ
lt_subplot(3,2,5); hold on;
title('learning vs. abs pos from targ');
X2=abs(X);
Y2=Y;
lt_regress(Y2, X2, 1, 0, 1, 1, 'b');
xlim([0 Xlim(2)+1]); lt_plot_zeroline;





% ===== diff type
lt_figure; hold on;
% 1) summary
lt_subplot(3,1,1); hold on;
title('diff type'); xlabel('pos rel targ'); ylabel('gener');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All);
color='r';

X=AllSylsStruct.PosRelTarg_All(inds); % pos rel targ
Y=AllSylsStruct.LearningRelTarg_all(inds);
plot(X, Y, 'o', 'Color',color);
Xlim=xlim;
lt_plot_zeroline;
set(gca, 'XTick', [min(X):max(X)])
% -anova
p=anovan(Y, {X});
lt_plot_pvalue(p, 'anova', 1);
% - plot mean
distanceslist=unique(X);
distanceslist=distanceslist(~isnan(distanceslist));
for i=1:length(distanceslist);
    dist=distanceslist(i);
    
    ymean=mean(Y(X==dist));
    ysem=lt_sem(Y(X==dist));
    lt_plot_bar(dist, ymean, {'Errors', ysem,'Color',color});
end
% 2) linear regression (hoff style) (preceding)
lt_subplot(3,2,3); hold on;
xlabel('position'); ylabel('generalization');
title('preceeding targ');
X2=X(X<0);
Y2=Y(X<0);
lt_regress(Y2, X2, 1, 0, 1, 1, color);
xlim([Xlim(1)-1 0]); lt_plot_zeroline;
% 3) linear regression (hoff style) (following)
lt_subplot(3,2,4); hold on;
title('following targ');
X2=X(X>0);
Y2=Y(X>0);
lt_regress(Y2, X2, 1, 0, 1, 1, color);
xlim([0 Xlim(2)+1]); lt_plot_zeroline;
% 4) learning vs. absolute pos from targ
lt_subplot(3,2,5); hold on;
title('learning vs. abs pos from targ');
X2=abs(X);
Y2=Y;
lt_regress(Y2, X2, 1, 0, 1, 1, color);
xlim([0 Xlim(2)+1]); lt_plot_zeroline;




%% MOTIF CONTINEUD
% - motif alone
inds=AllSylsStruct.Target_all==0;
clear group_anova
group_anova{1}=~isnan(AllSylsStruct.PosRelTarg_All(inds)); % 1 = same motif
Y = AllSylsStruct.LearningRelTarg_all(inds);

GroupNames={'samemotif'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);
[c, m, h, nms]=multcompare(stats);

disp(' ---  motif alone has no effect');
pause


% - motif position (same motif)
inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All) & ~AllSylsStruct.Similar_all==1;
clear group_anova
group_anova{1}=AllSylsStruct.PosRelTarg_All(inds); % position in same motif
% group_anova{2}=AllSylsStruct.Similar_all(inds);
Y = AllSylsStruct.LearningRelTarg_all(inds);

GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames, 'model','interaction');
[c, m, h, nms]=multcompare(stats);
disp('position (for same-type; diff type has no effect) has an effect by ANOVA, but post hoc nothing comes out');
pause



% === ANOVA
% type x motif(same?)

inds=AllSylsStruct.Target_all==0;
clear group_anova;
group_anova{1}=AllSylsStruct.Similar_all(inds);
group_anova{2}=~isnan(AllSylsStruct.PosRelTarg_All(inds)); % 1 = same motif
Y = AllSylsStruct.LearningRelTarg_all(inds);

GroupNames={'sametype', 'samemotif'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames, 'model','interaction');

[c, m, h, nms]=multcompare(stats, 'dimension', [1 2]);
disp(' ');
disp('NOTE: this suggests that only same-type and diff-motif is significant generalization');
pause;



% ==== why is generalization driven by different motif? (ordinal number?)
% correlate position 
% -- diff motif
inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.PosInMotif_regexp(inds)-AllSylsStruct.TargPosInMotif(inds);
Xname='PosInMot(rel targ)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
disp(' -- diff motif, no clear relationship between ordinal numb and gen');


% -- diff motif (abs value of distance from targ)
inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.PosInMotif_regexp(inds)-AllSylsStruct.TargPosInMotif(inds);
X=abs(X);
Xname='PosInMot(rel targ)(abs val)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
disp(' -- diff motif, no clear relationship between ordinal numb and gen');

% -- diff motif (abs value of distance from targ)
lt_figure; hold on;
% - same syl
lt_subplot(2,2,1); hold on;
title('same');
ylabel('gen'); xlabel('ordin pos rel targ (diff motif only');
inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All) & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.PosInMotif_regexp(inds)-AllSylsStruct.TargPosInMotif(inds);
X=abs(X);
Xname='PosInMot(rel targ)(abs val)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, 'b');
% -- diff syl
lt_subplot(2,2,2); hold on;
title('diff');
inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All) & AllSylsStruct.Similar_all==0;

X=AllSylsStruct.PosInMotif_regexp(inds)-AllSylsStruct.TargPosInMotif(inds);
X=abs(X);
Xname='PosInMot(rel targ)(abs val)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, 'r');



% ====
inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.PosInMotif_regexp(inds)-AllSylsStruct.TargPosInMotif(inds);
Xname='PosInMot(rel targ)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.Similar_all(inds);
Zname='type';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname);
disp(' -- same motif, no clear relationship between ordinal numb and gen');
pause


% ========== IS it because for same sequence (presim) on other motif?
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
clear group_anova;
group_anova{1}=AllSylsStruct.PreSimilar_all(inds);
group_anova{2}=~isnan(AllSylsStruct.PosRelTarg_All(inds)); % 1 = same motif
Y = AllSylsStruct.LearningRelTarg_all(inds);

GroupNames={'presim', 'samemotif'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames, 'model','interaction');

[c, m, h, nms]=multcompare(stats, 'dimension', [1 2]);
pause;


% ===== for all same-types, count number of presims on other motif vs. on
% same motif
sum(AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & AllSylsStruct.PreSimilar_all==1 & ~isnan(AllSylsStruct.PosRelTarg_All)) % presim + same motif

sum(AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & AllSylsStruct.PreSimilar_all==1 & isnan(AllSylsStruct.PosRelTarg_All)) % presim + diff motif

sum(AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & AllSylsStruct.PreSimilar_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All)) % prediff + same motif

sum(AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1 & AllSylsStruct.PreSimilar_all==0 & isnan(AllSylsStruct.PosRelTarg_All)) % prediff + diff motif



% =========
disp(' === CONCLUSION: generalization is greater for same-type vs diff type only for diff motif (for same motif is trend, but not significant). That effect');
disp('is driven by fact that there are proportionally more presim on diff motif (8/11 same types) vs. same motif (2/20 same types).');
disp(' === NO effect found for ordinal number on other motif. note that have not confirmed that the ordinal numbers I used for regexp motifs are ideal for this analysis');
disp(' ===  MOTIF posituion on same motif has effect for same-types, not diff types');


%% +++++++++++++ a bunch of stuff [not enough data]

inds=AllSylsStruct.Target_all==0;
clear group_anova;
group_anova{1}=AllSylsStruct.PreSimilar_all(inds); % presim
group_anova{2}=AllSylsStruct.Similar_all(inds); % type
% group_anova{3}=AllSylsStruct.BirdNum_all(inds);  % bird code

% tmp=AllSylsStruct.PosRelTarg_All(inds);
% tmp(isnan(tmp))=100; % recode diff motif as 100
% group_anova{4}=tmp; % motif position (or 100 for diff motif)
% 
Y = AllSylsStruct.LearningRelTarg_all(inds);

% GroupNames={'presim', 'syltype', 'birdnum', 'motifpos'};
GroupNames={'presim', 'syltype'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames, 'model','interaction');

[c, m, h, nms]=multcompare(stats, 'dimension', [1 2]);
pause;



%% +++++++++++++++++++++++++++++++++++++++++++++++ relationship between correlations, learning, and sequential position

% ===== correlation between sequential position and correlations?
inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All) & AllSylsStruct.Similar_all;


X=AllSylsStruct.Corr_motif_all(inds);
Xname='corr (motif)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.PosRelTarg_All(inds);
Zname='positionRelTarg?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');

% --- compare means
[~,~, stats]=anovan(Y, {Zgroups});
multcompare(stats)

disp('want to compare marginal means, but can"t because some values of motif position dont have enough data');
pause

% ===== correlation between sequential position (collapse into +1 and others) and correlations?
inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All) & AllSylsStruct.Similar_all;


X=AllSylsStruct.Corr_motif_all(inds);
Xname='corr (motif)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.PosRelTarg_All(inds);
Zgroups(Zgroups~=1)=100;
Zname='positionRelTarg?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'parallel lines');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');


disp('same-type, same-motif, after controlling for pos +1 versus rest, correlation between corr (motif) and learning is not significnat, and trending in opposite direction actually');
disp('in same analysis, effect of pos (+1 vs rest) is significant, suggesting that for same motif the effect of correlation vs. learning was driven solely by the fact that pos +1 syls learn more, and they happen to correlate more');
disp('DIFF type: positive correlation between corr and learning still holds even after taking into account pos +1 syl')
disp('in fact, pos+1 syl does not look any different from other syls');
pause


% ===== correlation between sequential position (collapse into +1 and others) and correlations?
inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All) & AllSylsStruct.Similar_all;


X=AllSylsStruct.Corr_motif_all(inds);
Xname='corr (motif)';
Y=AllSylsStruct.LearningRelTarg_all(inds);
Yname='gener';

Zgroups=AllSylsStruct.PosRelTarg_All(inds);
% Zgroups(Zgroups~=1)=100;
Zname='positionRelTarg?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'parallel lines');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');


%% =========== LMAN - does corr vs. song increase when LMAN inactivated?

% ===== 1) what happens to corr on average? (due to MUSC)
lt_figure; hold on;
numrows=3;

% -- same
lt_subplot(numrows,2,1); hold on;
title('same');
xlabel('corr-song, PBS');
ylabel('corr-song, MUSC');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
lt_plot_45degScatter(X, Y, 'b');

% -- diff
lt_subplot(numrows,2,2); hold on;
title('diff');
xlabel('corr-song, PBS');
ylabel('corr-song, MUSC');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
sylnames=AllSylsStruct.SylSingle_all(inds);
lt_plot_45degScatter(X, Y, 'r');
for i=1:length(sylnames);
    lt_plot_text(X(i)+0.05, Y(i), sylnames(i), 'b', 12);
end



% ---- diff effect depending on  motif [same motif]
lt_subplot(numrows,2,3); hold on;
title('diff [same motif]');
xlabel('corr-song, PBS');
ylabel('corr-song, MUSC');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
lt_plot_45degScatter(X, Y, 'r');

% ---- diff effect depending on  motif [diff motif]
lt_subplot(numrows,2,4); hold on;
title('diff [diff motif]');
xlabel('corr-song, PBS');
ylabel('corr-song, MUSC');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
lt_plot_45degScatter(X, Y, 'r');


% -- 
lt_subplot(numrows,2,5); hold on;
title('diff (handlab sim)');
xlabel('corr-song, PBS');
ylabel('corr-song, MUSC');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & strcmp(AllSylsStruct.SylSingle_all, 'b');

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
sylnames=AllSylsStruct.Syl_all(inds);
lt_plot_45degScatter(X, Y, 'r');
for i=1:length(sylnames);
    lt_plot_text(X(i)+0.05, Y(i), sylnames(i), 'b', 12);
end


% -- 
lt_subplot(numrows,2,6); hold on;
title('diff (handlab diff)');
xlabel('corr-song, PBS');
ylabel('corr-song, MUSC');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & ~strcmp(AllSylsStruct.SylSingle_all, 'b');

X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
sylnames=AllSylsStruct.Syl_all(inds);
lt_plot_45degScatter(X, Y, 'r');
for i=1:length(sylnames);
    lt_plot_text(X(i)+0.05, Y(i), sylnames(i), 'b', 12);
end

%% ===== does effect of LMAN on baseline correlation depend on some features?

lt_figure; hold on;
numrows=3;

% 1) ---- vs acoustic dist [all nontarg]
lt_subplot(numrows, 2, 1); hold on;
xlabel('acoust dist');
ylabel('MUSC corr - PBS corr');
title('all nontarg');

% - allnontarg
inds=AllSylsStruct.Target_all==0;
color='k';

X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline


% 1) ---- vs acoustic dist
lt_subplot(numrows, 2, 2); hold on;
xlabel('acoust dist');
ylabel('MUSC corr - PBS corr');
title('all nontarg');

% - same
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
color='b';

X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
lt_plot(X, Y, {'Color', color});
% - diff [hand sim]
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & strcmp(AllSylsStruct.SylSingle_all, 'b');

X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
plot(X, Y, 'or');
% - diff [hand diff]
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & ~strcmp(AllSylsStruct.SylSingle_all, 'b');

X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
lt_plot(X, Y, {'Color','r'});
% -
lt_plot_zeroline


% 2) ---- distributions
% - same & diff
lt_subplot(numrows, 2, 3); hold on;
xlabel('MUSC corr - PBS corr');
title('same & diff');

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
color='b';
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
YY=Y(~isnan(Y));
[~, X]=lt_plot_histogram(YY, '', 1, 0, '', 1, color);
line([0 0], ylim);

inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;
color='r';
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
YY=Y(~isnan(Y));
lt_plot_histogram(YY, X, 1, 0, '', 1, color);
line([0 0], ylim);


% - diff [hand same]
lt_subplot(numrows, 2, 4); hold on;
xlabel('MUSC corr - PBS corr');
title('diff(hand same)');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & strcmp(AllSylsStruct.SylSingle_all, 'b');
color='r';

Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
YY=Y(~isnan(Y));
lt_plot_histogram(YY, '', 1, 0, '', 1, color);
line([0 0], ylim);
% - diff [hand diff]
lt_subplot(numrows, 2, 5); hold on;
xlabel('MUSC corr - PBS corr');
title('diff(hand diff)');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0 & ~strcmp(AllSylsStruct.SylSingle_all, 'b');
color='r';

Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds);
YY=Y(~isnan(Y));
lt_plot_histogram(YY, '', 1, 0, '', 1, color);
line([0 0], ylim);


%% ------- COVARY WITH ACOUSTIC AND MOTIF [all nontarg] [continued from above]
inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds); % diff
Yname='corr (MUSC) - corr (PBS) [song]';
Zgroups=~isnan(AllSylsStruct.PosRelTarg_All(inds));
% Zgroups(Zgroups~=1)=100;
Zname='same motif?';


[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'parallel lines');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');


% ------- COVARY WITH ACOUSTIC AND hand lab same [diff]
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds)-AllSylsStruct.Corr_song_all(inds); % diff
Yname='corr (MUSC) - corr (PBS) [song]';
Zgroups=strcmp(AllSylsStruct.SylSingle_all(inds), 'b');
Zname='hand lab sim?';


[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'parallel lines');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');



disp('CONCLUSION: what moitf you are on does not affect how MUSC affects your corr with target')
disp('i.e. for both motifs average is increase in corr');
disp('CONCLUSION: for diff type, hand lab sim or diff does not affect the relationship between acoustic dist and effect of MUSC on corr');
disp('CONCLUSION: there is correlation (p=0.055) between effect of MUSC on corr and acoustic distance');



%% === relationship between corr and acoustic dist? (and relation to gener)

% 1) PBS (song)
inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.Corr_song_all(inds);
Yname='corr (song)';
Zgroups=AllSylsStruct.Similar_all(inds)==1;
Zname='same?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');


% 2) MUSC (song)
inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.LMAN.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
Yname='corr (song)';
Zgroups=AllSylsStruct.Similar_all(inds)==1;
Zname='same?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');



% 3) PBS (motif)
inds=AllSylsStruct.Target_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.Corr_motifsubsong_all(inds);
Yname='corr (motif)';
Zgroups=AllSylsStruct.Similar_all(inds)==1;
Zname='same?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');

% 3) PBS (diff)
inds=AllSylsStruct.Target_all==0 & isnan(AllSylsStruct.PosRelTarg_All);

X=AllSylsStruct.AcousticDist_all(inds);
Xname='acoust dist';
Y=AllSylsStruct.Corr_song_all(inds);
Yname='corr (song)';
Zgroups=AllSylsStruct.Similar_all(inds)==1;
Zname='same?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on');
[c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');



% ============== relation to gen (3d plot)
lt_figure; hold on;

% lt_subplot(2,2,1); hold on;
xlabel('acoustic dist');
ylabel('corr(song)');
zlabel('generalization');

% -- same
color='b';
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;

X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.Corr_song_all(inds);
Z=AllSylsStruct.LearningRelTarg_all(inds);
h=lt_plot_stem3(X, Y, Z, color, 0);


% -- diff
color='r';
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;

X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.Corr_song_all(inds);
Z=AllSylsStruct.LearningRelTarg_all(inds);
h=lt_plot_stem3(X, Y, Z, color, 1);


% lt_figure; hold on;
% 
% % lt_subplot(2,2,1); hold on;
% xlabel('acoustic dist');
% ylabel('corr(song)');
% zlabel('generalization');
% 
% % -- same
% color='b';
% inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
% 
% X=AllSylsStruct.AcousticDist_all(inds);
% Y=AllSylsStruct.Corr_song_all(inds);
% Z=AllSylsStruct.LearningRelTarg_all(inds);
% % -- plot pos and neg differently
% indstmp=Z>0;
% h=stem3(X(indstmp), Y(indstmp), Z(indstmp));
% set(h, 'Color',color, 'MarkerFaceColor', color)
% indstmp=Z<=0;
% h=stem3(X(indstmp), Y(indstmp), Z(indstmp));
% set(h, 'Color',color)
% 
% 
% % -- diff
% color='r';
% inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;
% 
% X=AllSylsStruct.AcousticDist_all(inds);
% Y=AllSylsStruct.Corr_song_all(inds);
% Z=AllSylsStruct.LearningRelTarg_all(inds);
% 
% stem3(X, Y, Z, 'Color',color);
% % -- plot pos and neg differently
% indstmp=Z>0;
% stem3(X(indstmp), Y(indstmp), Z(indstmp), 'Color',color, 'MarkerFaceColor', color);
% indstmp=Z<=0;
% stem3(X(indstmp), Y(indstmp), Z(indstmp), 'Color',color);
% 
% x=xlim; x=linspace(x(1), x(2), 10);
% y=ylim; y=linspace(y(1), y(2), 10);
% z=zeros(length(x), length(y));
% mesh(x, y, z);


%% ====== MUSC corrs better predict generalization than PBS corrs?
lt_figure; hold on;

% 1) all nontarg
lt_subplot(3,2,1); hold on;
xlabel('corr(song)');
ylabel('gener');
lt_plot_annotation(1, 'bk: PBS; rd: MUSC [all nontarg]', 'k');
inds=AllSylsStruct.Target_all==0;
% - PBS
color='k';
X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);

% 1) all nontarg
lt_subplot(3,2,2); hold on;
xlabel('corr(motif)');
ylabel('gener');
lt_plot_annotation(1, 'bk: PBS; rd: MUSC [all nontarg, same motif]', 'k');
inds=AllSylsStruct.Target_all==0 & ~isnan(AllSylsStruct.PosRelTarg_All);
% - PBS
color='k';
X=AllSylsStruct.Corr_motif_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.Corr_motif_all_MUSC(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);


% 1) diff type
lt_subplot(3,2,3); hold on;
xlabel('corr(song)');
ylabel('gener');
title('bk: PBS; rd: MUSC [diff type]');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;
% - PBS
color='k';
X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);


% 1) same type
lt_subplot(3,2,4); hold on;
xlabel('corr(song)');
ylabel('gener');
title('bk: PBS; rd: MUSC [same type]');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
% - PBS
color='k';
X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);

% 1) same type [hand lab]
lt_subplot(3,2,5); hold on;
xlabel('corr(song)');
ylabel('gener');
title('bk: PBS; rd: MUSC [st; hand lab]');
inds=AllSylsStruct.Target_all==0 & strcmp(AllSylsStruct.SylSingle_all, 'b');
% - PBS
color='k';
X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);

% 1) diff type [hand lab]
lt_subplot(3,2,6); hold on;
xlabel('corr(song)');
ylabel('gener');
title('bk: PBS; rd: MUSC [diff; hand lab]');
inds=AllSylsStruct.Target_all==0 & ~strcmp(AllSylsStruct.SylSingle_all, 'b');
% - PBS
color='k';
X=AllSylsStruct.Corr_song_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.Corr_song_all_MUSC(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);


disp('CONCLUSION: LMAN inactivation did not increase the ability of corr to predict generalization - did slightly but almost certain is not significant');
disp('that is true for all nontarg, for same motif or all motif, for diff types');

%% ===== INACTIVATION - RELATIONSHIP BETWEEN ACOUSTIC AND LEARNING CHANGES?
lt_figure; hold on;

% 1) sametype
lt_subplot(3,2,1); hold on;
xlabel('acoustic dist');
ylabel('gener');
lt_plot_annotation(1, 'bk: PBS; rd: MUSC [same]', 'k');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
% - PBS
color='k';
X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.AcousticDist_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);

% 1) diff
lt_subplot(3,2,2); hold on;
xlabel('acoustic dist');
ylabel('gener');
lt_plot_annotation(1, 'bk: PBS; rd: MUSC [diff]', 'k');
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;
% - PBS
color='k';
X=AllSylsStruct.AcousticDist_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);
% - MUSC
color='r';
X=AllSylsStruct.LMAN.AcousticDist_all(inds);
Y=AllSylsStruct.LearningRelTarg_all(inds);
lt_regress(Y, X, 1, 0, 1, 1, color);


%% === plot distributions of acoustic distances (LMAN active and ianctive)

% 1) PBS
lt_subplot(3,2,3); hold on;
xlabel('acoustic dist');
title('PBS');
Xcenters=0:0.2:8;
% -- same
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
color='b';

Y=AllSylsStruct.AcousticDist_all(inds);
lt_plot_histogram(Y, Xcenters, 1, 0, '', 1, color);
% -- diff
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;
color='r';

Y=AllSylsStruct.AcousticDist_all(inds);
lt_plot_histogram(Y, Xcenters, 1, 0, '', 1, color);

% 1) PBS
lt_subplot(3,2,4); hold on;
xlabel('acoustic dist');
title('MUSC');
Xcenters=0:0.2:8;
% -- same
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==1;
color='b';

Y=AllSylsStruct.LMAN.AcousticDist_all(inds);
lt_plot_histogram(Y, Xcenters, 1, 0, '', 1, color);
% -- diff
inds=AllSylsStruct.Target_all==0 & AllSylsStruct.Similar_all==0;
color='r';

Y=AllSylsStruct.LMAN.AcousticDist_all(inds);
lt_plot_histogram(Y, Xcenters, 1, 0, '', 1, color);

