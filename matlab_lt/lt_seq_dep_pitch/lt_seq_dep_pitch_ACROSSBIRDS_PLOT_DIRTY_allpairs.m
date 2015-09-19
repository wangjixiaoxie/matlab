function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_allpairs(SeqDepPitch_AcrossBirds, PARAMS)

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% %% [DIRTY!!!!!!!!!!!!] [PAIRWISE STATISTICS] [PLOT ACOUSTIC DISTANCE REL TO DISTANCE IN MOTIF (NOTE: THIS USES THE LEARNING CODE ABOVE, AND IS NOT GREAT. LITERALLY PUTS ACOUSTIC DISTANCE IN THE LAERNING VARIABLE
clear syl;
write_corr=0;
write_syl=0;

% ONLY ONE OF THESE CAN BE 1
plot_acoustic=0;
plot_corr_song=1;
plot_corr_motif=0;


lt_figure; hold on;

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
if plot_corr_song==1;
ylim([-1 1]);
ylabel('corr (song)');
end

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
if plot_corr_song==1;
ylim([-1 1]);
ylabel('corr (song)');
end
xlabel('motif distance');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
if plot_corr_song==1;
ylim([-1 1]);
ylabel('corr (song)');
end
xlabel('motif distance');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
if plot_corr_song==1;
ylim([-1 1]);
ylabel('corr (song)');
end

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
xlabel('motif distance');
if plot_corr_song==1;
ylim([-1 1]);
ylabel('corr (song)');
end

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
xlabel('motif distance');
if plot_corr_song==1;
ylim([-1 1]);
ylabel('corr (song)');
end

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        %         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        for j=1:length(syllist);
            syl1=syllist{j};
            % skip is last syl
            if j==length(syllist);
                continue;
            end
            
            for jj=j+1:length(syllist);
                syl2=syllist{jj};
                
                % skip if this is a target
                %             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                %                 continue;
                %             end
                
                % ==== collect data
                % acoustic distance
                fv1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl1).fv_baseline_zscore_mean;
                fv2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_zscore_mean;
                acoustdist=sqrt(sum((fv1-fv2).^2));

                % correlations
                try
                    corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                catch err
                    corr_motif=nan;
                end
                corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                
                
                % put into variable that will plot
                if plot_acoustic==1;
                    learning_rel_targ=acoustdist;
                elseif plot_corr_song==1;
                    learning_rel_targ=corr_song;
                elseif plot_corr_motif==1;
                    learning_rel_targ=corr_motif;
                end
                Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
                
                
                % - motif distance
                motifnum1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_motifnum;
                motifnum2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_motifnum;
                if motifnum1~=motifnum2;
                    motif_distance=nan;
                else
                    motifpos1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_PosInMotif_thissyl;
                    motifpos2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_PosInMotif_thissyl;
                    
                    motif_distance=abs(motifpos1-motifpos2);
                end
                
                MotifDist_AllExpts=[MotifDist_AllExpts motif_distance];
                
                %             % acoustic distance
                %             eucl_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
                
                % is similar?
                lowersyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).single_syl;
                lowersyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).single_syl;
                if lowersyl1==lowersyl2;
                    is_similar=1;
                else
                    is_similar=0;
                end
                
                Similar_AllExpts=[Similar_AllExpts is_similar];
                
                
                % predecing syllable similar?
                presyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).preceding_syl;
                if length(presyl1)>1;
                    presyl1=regexp(presyl1, '[A-Z]', 'match');
                    presyl1=lower(presyl1);
                end
                presyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).preceding_syl;
                if length(presyl2)>1;
                    presyl2=regexp(presyl2, '[A-Z]', 'match');
                    presyl2=lower(presyl2);
                end
                
                if strcmp(presyl1, presyl2);
                    preceding_similar=1;
                else
                    preceding_similar=0;
                end
                PreSylSimilar_AllExpts=[PreSylSimilar_AllExpts preceding_similar];
                
                
                
                % subplots, for those on same and those on different motif
                if isnan(motif_distance); % then not on same motif
                    subplot(2,4,1);
                    
                    %             ==== plot
                    X=-0.5+rand;
                    if is_similar==1;
                        lt_plot(X, learning_rel_targ, {'Color','b'});
                    else
                        lt_plot(X, learning_rel_targ, {'Color','r'});
                    end
                    
                    % if presyl is similar, outline in black
                    if preceding_similar==1;
                        plot(X, learning_rel_targ, 'ok','MarkerSize',7, 'LineWidth', 2.5);
                    end
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                else
                    % == plot in "in target" subplot
                    subplot(2,4,2);
                    %     --- plot
                    X=-0.5+rand;
                    if is_similar==1;
                        lt_plot(X, learning_rel_targ, {'Color','b'});
                    else
                        lt_plot(X, learning_rel_targ, {'Color','r'});
                    end
                    
                    % if presyl is similar, outline in black
                    if preceding_similar==1;
                        plot(X, learning_rel_targ, 'ok','MarkerSize',7,'LineWidth', 2.5);
                    end
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                    
                    
                    % === Plot based on motif distance
                    subplot(2,4,3:4);
                    %  -- plot
                    X=motif_distance-0.2+0.4*rand;
                    if is_similar==1;
                        lt_plot(X, learning_rel_targ, {'Color','b'});
                    else
                        lt_plot(X, learning_rel_targ, {'Color','r'});
                    end
                    
                    % if presyl is similar, outline in black
                    if preceding_similar==1;
                        plot(X, learning_rel_targ, 'ok','MarkerSize',7, 'LineWidth', 2.5);
                    end
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                    
                end
            end
        end
    end
end

if write_syl==1;
    disp('bird syl1 syl2')
end



% ==== PLOT ACROSS EXPERIMENTS
% ---- Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.SimilarToTarg=similar;

% all
learning_mean=nanmean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);


% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.similar.LearningRelTarg=Y;

% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);


% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=nanmean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.different.LearningRelTarg=Y;
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=nanmean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);




% ---- Same Motifs (all combined);
subplot(2,4,6); hold on;
inds=~isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

% all
learning_mean=mean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);

% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);



% ---- Same Motifs (based on distance
subplot(2,4,7:8); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);
preceding=PreSylSimilar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SimilarToTarg=similar;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.DistanceInMotif=motif_distance;


distances_that_exist=unique(motif_distance);
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    Y_preceding=preceding(inds2);
    Y_learning=learning(inds2);
    Y_similar=similar(inds2);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.distances_that_exist=distances_that_exist;
    
    % -- all
    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.all.LearningRelTarg=Y_learning;
    
    % -- similar
    inds3=Y_similar==1;
    Y_pre=Y_preceding(inds3);

    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.LearningRelTarg=Y_learning(inds3);
    
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerSize', 8);
    
    
    % -- diff
    inds3=Y_similar==0;
    Y_pre=Y_preceding(inds3);
    
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));

    errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.LearningRelTarg=Y_learning(inds3);
        
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerSize', 8);

    
    
    
    if plot_acoustic==1;
        lt_subtitle('ALL PAIRS OF SYLS - ACOUSTIC vs. MOTIF POSITION')
    elseif plot_corr_song==1;
        lt_subtitle('ALL PAIRS OF SYLS - CORR (song) vs. MOTIF POSITION')
    elseif plot_corr_motif==1;
        lt_subtitle('ALL PAIRS OF SYLS - COPR (motif) vs. MOTIF POSITION')
    end
    
end

%% %% [REPEATS ABOVE, IGNORE (DO CORR AND ACOUSTIC ONCE EACH) [DIRTY!!!!!!!!!!!!] [PAIRWISE STATISTICS] [PLOT CORR DISTANCE REL TO DISTANCE IN MOTIF (NOTE: THIS USES THE LEARNING CODE ABOVE, AND IS NOT GREAT. LITERALLY PUTS ACOUSTIC DISTANCE IN THE LAERNING VARIABLE
clear syl;
write_corr=0;
write_syl=0;

% ONLY ONE OF THESE CAN BE 1
plot_acoustic=0;
plot_corr_song=1;
plot_corr_motif=0;


lt_figure; hold on;

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-1 1]);
ylabel('learning (rel to target)');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
% ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
ylabel('learning (rel to target)');

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        %         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        for j=1:length(syllist);
            syl1=syllist{j};
            % skip is last syl
            if j==length(syllist);
                continue;
            end
            
            for jj=j+1:length(syllist);
                syl2=syllist{jj};
                
                % skip if this is a target
                %             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                %                 continue;
                %             end
                
                % ==== collect data
                % acoustic distance
                fv1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl1).fv_baseline_zscore_mean;
                fv2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_zscore_mean;
                acoustdist=sqrt(sum((fv1-fv2).^2));

                % correlations
                try
                    corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                catch err
                    corr_motif=nan;
                end
                corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                
                
                % put into variable that will plot
                if plot_acoustic==1;
                    learning_rel_targ=acoustdist;
                elseif plot_corr_song==1;
                    learning_rel_targ=corr_song;
                elseif plot_corr_motif==1;
                    learning_rel_targ=corr_motif;
                end
                Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
                
                
                % - motif distance
                motifnum1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_motifnum;
                motifnum2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_motifnum;
                if motifnum1~=motifnum2;
                    motif_distance=nan;
                else
                    motifpos1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_PosInMotif_thissyl;
                    motifpos2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_PosInMotif_thissyl;
                    
                    motif_distance=abs(motifpos1-motifpos2);
                end
                
                MotifDist_AllExpts=[MotifDist_AllExpts motif_distance];
                
                %             % acoustic distance
                %             eucl_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
                
                % is similar?
                lowersyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).single_syl;
                lowersyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).single_syl;
                if lowersyl1==lowersyl2;
                    is_similar=1;
                else
                    is_similar=0;
                end
                
                Similar_AllExpts=[Similar_AllExpts is_similar];
                
                
                % predecing syllable similar?
                presyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).preceding_syl;
                if length(presyl1)>1;
                    presyl1=regexp(presyl1, '[A-Z]', 'match');
                    presyl1=lower(presyl1);
                end
                presyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).preceding_syl;
                if length(presyl2)>1;
                    presyl2=regexp(presyl2, '[A-Z]', 'match');
                    presyl2=lower(presyl2);
                end
                
                if strcmp(presyl1, presyl2);
                    preceding_similar=1;
                else
                    preceding_similar=0;
                end
                PreSylSimilar_AllExpts=[PreSylSimilar_AllExpts preceding_similar];
                
                
                
                % subplots, for those on same and those on different motif
                if isnan(motif_distance); % then not on same motif
                    subplot(2,4,1);
                    
                    %             ==== plot
                    X=-0.5+rand;
                    if is_similar==1;
                        lt_plot(X, learning_rel_targ, {'Color','b'});
                    else
                        lt_plot(X, learning_rel_targ, {'Color','r'});
                    end
                    
                    % if presyl is similar, outline in black
                    if preceding_similar==1;
                        plot(X, learning_rel_targ, 'ok','MarkerSize',7, 'LineWidth', 2.5);
                    end
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                else
                    % == plot in "in target" subplot
                    subplot(2,4,2);
                    %     --- plot
                    X=-0.5+rand;
                    if is_similar==1;
                        lt_plot(X, learning_rel_targ, {'Color','b'});
                    else
                        lt_plot(X, learning_rel_targ, {'Color','r'});
                    end
                    
                    % if presyl is similar, outline in black
                    if preceding_similar==1;
                        plot(X, learning_rel_targ, 'ok','MarkerSize',7,'LineWidth', 2.5);
                    end
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                    
                    
                    % === Plot based on motif distance
                    subplot(2,4,3:4);
                    %  -- plot
                    X=motif_distance-0.2+0.4*rand;
                    if is_similar==1;
                        lt_plot(X, learning_rel_targ, {'Color','b'});
                    else
                        lt_plot(X, learning_rel_targ, {'Color','r'});
                    end
                    
                    % if presyl is similar, outline in black
                    if preceding_similar==1;
                        plot(X, learning_rel_targ, 'ok','MarkerSize',7, 'LineWidth', 2.5);
                    end
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                    
                end
            end
        end
    end
end

if write_syl==1;
    disp('bird syl1 syl2')
end



% ==== PLOT ACROSS EXPERIMENTS
% ---- Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.SimilarToTarg=similar;

% all
learning_mean=nanmean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);


% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.similar.LearningRelTarg=Y;

% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);


% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=nanmean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.different.LearningRelTarg=Y;
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=nanmean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);




% ---- Same Motifs (all combined);
subplot(2,4,6); hold on;
inds=~isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

% all
learning_mean=mean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);

% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);



% ---- Same Motifs (based on distance
subplot(2,4,7:8); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);
preceding=PreSylSimilar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SimilarToTarg=similar;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.DistanceInMotif=motif_distance;


distances_that_exist=unique(motif_distance);
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    Y_preceding=preceding(inds2);
    Y_learning=learning(inds2);
    Y_similar=similar(inds2);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.distances_that_exist=distances_that_exist;
    
    % -- all
    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.all.LearningRelTarg=Y_learning;
    
    % -- similar
    inds3=Y_similar==1;
    Y_pre=Y_preceding(inds3);

    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.LearningRelTarg=Y_learning(inds3);
    
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerSize', 8);
    
    
    % -- diff
    inds3=Y_similar==0;
    Y_pre=Y_preceding(inds3);
    
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));

    errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.LearningRelTarg=Y_learning(inds3);
        
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerSize', 8);

    
    
    
    if plot_acoustic==1;
        lt_subtitle('ALL PAIRS OF SYLS - ACOUSTIC vs. MOTIF POSITION')
    elseif plot_corr_song==1;
        lt_subtitle('ALL PAIRS OF SYLS - CORR (song) vs. MOTIF POSITION')
    elseif plot_corr_motif==1;
        lt_subtitle('ALL PAIRS OF SYLS - COPR (motif) vs. MOTIF POSITION')
    end
    
end

