function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbase(SeqDepPitch_AcrossBirds, PARAMS, same_type_thr, NullControl_splitday)
%% analysis of effect of muscimol on baseline pitch

% NullControl_splitday=1; % then isntead of taking LMAN days takes any experiment with PBS days with data in the afternoon, and splits data into late (>2pm, MUSC) and early (<1pm, PBS), to see if below effeects
% are clearly due to time of day. DRAWBACK - not using same experiments and
% only some birds are the same


%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
AcousVec_PBS=[];
Pitch_MUSC_minus_PBS=[];
Pitch_PBSraw=[];
PitchSTD_PBS=[];
ExptCounter=[];
Birdnames={};
Exptnames={};
Syllables={};

% === pairs that cannot be obtained in next section( i.e. must be gotten
% from raw data);
CorrSong_PBS_pairs=[];
CorrSong_MUSC_pairs=[];
CorrMotif_PBS_pairs=[];
CorrMotif_MUSC_pairs=[];
PairedAcousticDist=[];
PairedAcousticDist_MUSC=[];

exptcount=0;
counttmp=1;
if NullControl_splitday==0;
    % copy strcuture, save backup.
    filter = 'LMAN';
    [SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
    
    
    % ----- COLLECT data
    
    % acoustic vector (PBS)
    % acoustic vect (MUSC)
    % corr matrix (PBS)
    % corr (MUSC)
    % pitch shift caused by MUSC
    % baseline pitch during PBS
    % bird/expt ID
    % Syl, Bird, Expt
    
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexperiments;
            exptcount=exptcount+1;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            for j=1:length(SylsUnique)
                syl=SylsUnique{j};
                
                % --- extract syl stuff
                baseline_pbs=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                pitchstd_pbs=std(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                
                baseline_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                baseline_musc_minus_pbs=baseline_musc-baseline_pbs;
                
                acoustic_vec_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore_mean;
                acoustic_vec_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_zscore_mean;
                % ---- extract paired stuff
                for jj=j+1:length(SylsUnique);
                    syl2=SylsUnique{jj};
                    
                    % --- sanity check - so that matches next section
                    % paired code
                    acoustictmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_zscore_mean;
                    acousticdist=sqrt(sum((acoustic_vec_PBS-acoustictmp).^2));
                    
                    acoustictmp_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl2).fv_baseline_zscore_mean;
                    acousticdist_MUSC=sqrt(sum((acoustic_vec_MUSC-acoustictmp_musc).^2));
                    
                    corrsongPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                    corrsongMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                    
                    corrmotifPBS=nan;
                    corrmotifMUSC=nan;
                    
                    if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS, 'motif_by_motif');
                        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, syl2);
                            
                            corrmotifPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                            corrmotifMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                        end
                    end
                    if isempty(corrmotifPBS)
                        keyboard
                    end
                    % ===== OUTPUT
                    CorrMotif_PBS_pairs=[CorrMotif_PBS_pairs corrmotifPBS]; 
                    CorrMotif_MUSC_pairs=[CorrMotif_MUSC_pairs corrmotifMUSC];
                    CorrSong_PBS_pairs=[CorrSong_PBS_pairs corrsongPBS]; 
                    CorrSong_MUSC_pairs=[CorrSong_MUSC_pairs corrsongMUSC];
                    PairedAcousticDist=[PairedAcousticDist acousticdist];
                    counttmp=[counttmp counttmp(end)+1];
                    
                    PairedAcousticDist_MUSC=[PairedAcousticDist_MUSC acousticdist_MUSC];
                    
                end
                
                
                % ======= COLLECT
                AcousVec_PBS=[AcousVec_PBS acoustic_vec_PBS'];
                Pitch_MUSC_minus_PBS=[Pitch_MUSC_minus_PBS baseline_musc_minus_pbs];
                Pitch_PBSraw=[Pitch_PBSraw baseline_pbs];
                PitchSTD_PBS=[PitchSTD_PBS pitchstd_pbs];
                ExptCounter=[ExptCounter exptcount];
                Birdnames=[Birdnames birdname];
                Exptnames=[Exptnames exptname];
                Syllables=[Syllables syl];
                
            end
        end
    end
    
    
else
    
    figcount=1;
    subplotrows=5;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    
    
    
    NumBirds=length(SeqDepPitch_AcrossBirds.birds);
    % --- then get mock data
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexperiments;
            
            % === skip if is LMAN experiment
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                continue
            end
            
            
            
            exptcount=exptcount+1;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            for j=1:length(SylsUnique)
                syl=SylsUnique{j};
                
                % determine which tvals are mock "PBS" and which are mock
                % "MUSC"
                baseline_tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals;
                [~, tmp]=lt_convert_datenum_to_hour(baseline_tvals);
                baseline_tvals=tmp.hours;
                
                
                % -- plot distrubtion of hours
                if j==1
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title([birdname '-' exptname '-' syl]);
                    xlabel('time of baseline datapt');
                    lt_plot_histogram(baseline_tvals, '', 1, 0, 0, 1, 'k');
                    line([13 13], ylim, 'Color','r'); line([14 14], ylim); line([17.5 17.5], ylim);
                end
                
                
                % ---- mock inds
                inds_PBS=find(baseline_tvals<13);
                inds_MUSC=find(baseline_tvals>14 & baseline_tvals<17.5);
                
                % ===== extract syl stuff
                baseline_ffvals_all=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF;
                
                baseline_pbs=mean(baseline_ffvals_all(inds_PBS));
                pitchstd_pbs=std(baseline_ffvals_all(inds_PBS));
                
                baseline_musc=mean(baseline_ffvals_all(inds_MUSC));
                
                
                
                
                baseline_musc_minus_pbs=baseline_musc-baseline_pbs;
                
                acoustic_vec_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore_mean;
                
                
                
                
                % ======= COLLECT
                AcousVec_PBS=[AcousVec_PBS acoustic_vec_PBS'];
                Pitch_MUSC_minus_PBS=[Pitch_MUSC_minus_PBS baseline_musc_minus_pbs];
                Pitch_PBSraw=[Pitch_PBSraw baseline_pbs];
                PitchSTD_PBS=[PitchSTD_PBS pitchstd_pbs];
                ExptCounter=[ExptCounter exptcount];
                Birdnames=[Birdnames birdname];
                Exptnames=[Exptnames exptname];
                Syllables=[Syllables syl];
                
            end
        end
    end
    
    
end

% ===== other stuff
Pitch_MUSC_minus_PBS_zscore=Pitch_MUSC_minus_PBS./PitchSTD_PBS;

%% +++++++++++++ PLOTS

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];




%% ====== 1) distribution of shift caused by muscimol
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('distr of baseline musc shifts');

lt_plot_histogram(Pitch_MUSC_minus_PBS, '', 1, 1, 0, 1, 'k');
p=signrank(Pitch_MUSC_minus_PBS);
lt_plot_pvalue(p, 'signrank')

% --- use zscore
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('distr of baseline musc shifts');
ylabel('zscore of PBS across rends/days');
lt_plot_histogram(Pitch_MUSC_minus_PBS_zscore, '', 1, 1, 0, 1, 'k');
p=signrank(Pitch_MUSC_minus_PBS_zscore);
lt_plot_pvalue(p, 'signrank')


% --- does shift correlate with baseline std?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('musc shifts vs. baseline std (across rends and days)');

plot(PitchSTD_PBS, Pitch_MUSC_minus_PBS, 'ok');

% --- does shift correlate with baseline std?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('abs(musc shifts) vs. baseline std (across rends and days)');

plot(PitchSTD_PBS, abs(Pitch_MUSC_minus_PBS), 'ok');




%% ======= 2) if you are more similar structurally, do you have more similar musc shift?
%                 AcousVec_PBS=[AcousVec_PBS acoustic_vec_PBS'];
%                 Pitch_MUSC_minus_PBS=[Pitch_MUSC_minus_PBS baseline_musc_minus_pbs];
%                 Pitch_PBSraw=[Pitch_PBSraw baseline_pbs];
%                 PitchSTD_PBS=[PitchSTD_PBS pitchstd_pbs];
%                 ExptCounter=[ExptCounter exptcount];
%                 Birdnames=[Birdnames birdname];
%                 Exptnames=[Exptnames exptname];
%                 Syllables=[Syllables syl];

AcousticDist_pairs=[];
DiffInPitchShiftCausedByMusc_pairs=[];
similarDirOfShift_pairs=[];
changeInAbsSeparationDuringMUSC_pairs=[];
separationPBS_paired=[];
separationMUSC_paired=[];
PitchShiftBothSyls_paired={};


% --- for each experiment, get all pairwise 1) acoustic distance, and 2)
% similarity of shift
NumExpts=max(ExptCounter);

for i=1:NumExpts;
    Inds=find(ExptCounter==i);
    
    % ---- get all pairwise stuff
    for j=1:length(Inds);
        ind1=Inds(j);
        
        for jj=j+1:length(Inds);
            ind2=Inds(jj);
            
            acousticDist=sqrt(sum((AcousVec_PBS(:, ind1)-AcousVec_PBS(:, ind2)).^2));
            muscEffectDiff=abs(Pitch_MUSC_minus_PBS(ind1)-Pitch_MUSC_minus_PBS(ind2));
            similarDirOfShift=sign(Pitch_MUSC_minus_PBS(ind1))==sign(Pitch_MUSC_minus_PBS(ind2));
            
            separationPBS=Pitch_PBSraw(ind1)-Pitch_PBSraw(ind2);
            separationMUSC=(Pitch_PBSraw(ind1)+Pitch_MUSC_minus_PBS(ind1)) - (Pitch_PBSraw(ind2)+Pitch_MUSC_minus_PBS(ind2));
            changeInAbsSeparationDuringMUSC=abs(separationMUSC)-abs(separationPBS);
            
            
            % ==== collect pairwise stats
            PitchShiftBothSyls_paired=[PitchShiftBothSyls_paired [Pitch_MUSC_minus_PBS(ind1) Pitch_MUSC_minus_PBS(ind2)]];
            AcousticDist_pairs=[AcousticDist_pairs acousticDist];
            DiffInPitchShiftCausedByMusc_pairs=[DiffInPitchShiftCausedByMusc_pairs muscEffectDiff];
            similarDirOfShift_pairs=[similarDirOfShift_pairs similarDirOfShift];
            changeInAbsSeparationDuringMUSC_pairs=[changeInAbsSeparationDuringMUSC_pairs changeInAbsSeparationDuringMUSC];
            separationPBS_paired=[separationPBS_paired separationPBS];
            separationMUSC_paired=[separationMUSC_paired separationMUSC];
            
            
            %            disp([num2str(ind1) ' - '  num2str(ind2)])
            
        end
    end
end

% == confirm that things match from previuos code
if NullControl_splitday==0
    lt_figure; hold on;
    title('checking code match - must be 45deg');
    lt_plot_45degScatter(AcousticDist_pairs, PairedAcousticDist);
end

%% ====== PLOT


figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(AcousticDist_pairs>0);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(AcousticDist_pairs<same_type_thr);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
color='r';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(AcousticDist_pairs>=same_type_thr);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;



% === 1) diff in effect of musc vs. acoustic distance
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] diff in effect of musc vs. acoustic distance');
xlabel('acoustic dist')
ylabel('abs diff in effect of musc (hz)');

lt_regress(DiffInPitchShiftCausedByMusc_pairs, AcousticDist_pairs, 1, 0, 1, 1, 'k');


% === 2) similarity of sign of musc effect vs. acoustic distance
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] similarity of sign of musc effect vs. acoustic distance');
xlabel('acoustic dist')
ylabel('similarity of sign of musc effect');

lt_regress(similarDirOfShift_pairs, AcousticDist_pairs, 1, 0, 1, 1, 'k');


% === 2) similarity of sign of musc effect vs. acoustic distance
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] similarity of sign of musc effect vs. acoustic distance');
xlabel('acoustic dist')
ylabel('similarity of sign of musc effect');

[~, inds]=sort(AcousticDist_pairs);
X=AcousticDist_pairs(inds);
Y=similarDirOfShift_pairs(inds);

plot(X, Y, 'ok');
Yrun=lt_running_stats(Y, 50);
Xrun=lt_running_stats(X, 50);

shadedErrorBar(Xrun.Mean, Yrun.Mean, Yrun.SEM, {'Color','b'}, 1);


% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=AcousticDist_pairs<same_type_thr;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=AcousticDist_pairs>=same_type_thr;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

Y=changeInAbsSeparationDuringMUSC_pairs;

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=AcousticDist_pairs<same_type_thr;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'b');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds=AcousticDist_pairs<same_type_thr;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=AcousticDist_pairs>=same_type_thr;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'r');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds=AcousticDist_pairs>=same_type_thr;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'r');



% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] ');
xlabel('acoustic dist')
ylabel('change in separation of raw pitch during musc vs pbs');

lt_regress(changeInAbsSeparationDuringMUSC_pairs, AcousticDist_pairs, 1, 0, 1, 1, 'k');

%% ===== MUSC drive same-type to be more correlated? [song corr]
% --- 45 deg (same)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');
inds=AcousticDist_pairs<same_type_thr;
color='b';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);

% --
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');
inds=AcousticDist_pairs>=same_type_thr;
color='r';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);

% -----
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('corr(MUSC)-corr(PBS)');
X=-1:0.05:1;
% diff
inds=AcousticDist_pairs>=same_type_thr;
color='r';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);
% same
inds=AcousticDist_pairs<same_type_thr;
color='b';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

%% ===== MUSC drive same-type to be more correlated? [motif corr]
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



%% ==== effect of MUSC on corr correlate with effect of MUSC on separation?
% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('corr(MUSC) - corr(PBS)[song]');
ylabel('change in abs separation');
inds=AcousticDist_pairs<same_type_thr;
color='b';

X=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

% --- diff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('corr(MUSC) - corr(PBS)[song]');
ylabel('change in abs separation');
inds=AcousticDist_pairs>=same_type_thr;
color='r';

X=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)


disp('Effect of MUSC on pairwise corr did not predict effect of MUSC on separation [both sim and diff]');


%% === change in corr versus corr
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

%% ==== distributions of correlations (MUSC and PBS);

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

% --
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
inds=AcousticDist_pairs<same_type_thr;
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


% --
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
inds=AcousticDist_pairs>=same_type_thr;
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


%% ==== effect of musc on acoustic separation,

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
inds=find(AcousticDist_pairs>0);

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same]');
color='b';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
inds=find(AcousticDist_pairs<PARAMS.SylClassify.SylDistCutoff);

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff]');
color='r';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
inds=find(AcousticDist_pairs>PARAMS.SylClassify.SylDistCutoff);

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


%% ==== change in acoustic separation (vs. MUSC acoustic sep)
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
inds=find(AcousticDist_pairs>0);

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same]');
color='b';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
inds=find(AcousticDist_pairs<PARAMS.SylClassify.SylDistCutoff);

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff]');
color='r';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
inds=find(AcousticDist_pairs>PARAMS.SylClassify.SylDistCutoff);

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

