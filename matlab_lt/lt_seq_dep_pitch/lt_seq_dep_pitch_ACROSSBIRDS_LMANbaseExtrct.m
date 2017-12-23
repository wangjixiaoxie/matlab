function OUTSTRUCT=lt_seq_dep_pitch_ACROSSBIRDS_LMANbaseExtrct(SeqDepPitch_AcrossBirds, ...
    PARAMS, same_type_thr, NullControl_splitday, useHandLabForSametype, ...
    onlyUseFirstExpt, recalcValues, reCalcOldMethodKeepingNanSongs, UseMotifByMotifCorr)
%%

% ========= for recalculating song corr using regexp data (as origianlly
% done)
OnlyKeepDaysWithMUSC = 1; % default = 1 (i.e. matched days)
% reCalcOldMethodKeepingNanSongs = 1;
CorrSubtractDayMean=1; % default = 1; i.e. subtract day mean for all corr analyses.

%% lt 12/11/17 - separated code into data extraction (this) and separate plotting functions (others)

OUTSTRUCT = struct;

%% 12/10/17
UseLMANDatForSplit = 1; % i.e. use same dataset for both data and for negative control (split analsysi)
% if 1, uses only LMAN. if 0, uses only not LMAN.

% NOTE (IMPORTANT) - current progress for null control doing split - ran
% into issue in that LMAN experiments I did not label afetrnoon, even for
% the PBS days. Would like to use those days as control (afternooon vs.
% momring) but cannot. Solution is to look at other non-LMAN experiments
% for those birds? What other solutions? Look at raw plots of sample size
% versus hour. Need to implement changes at that lcoation.



%% lt 11/30/17 - use hand label instead of acousti threshodl to determine if similar
% NOTE: in some cases will do both hand lab and thresholded, in others will
% decide based on this variable.
% useHandLabForSametype=1;

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
SingleSylAll = {};

PitchTvalsPBS = {};
PitchTvalsMUSC = {};
PitchFFvalsPBS = {};
PitchFFvalsMUSC = {};

LearnHzAll = [];
LearnHzTargAll = [];
IsTargAll = [];
TargLearnDir = [];
DistFromTarg = [];

% === pairs that cannot be obtained in next section( i.e. must be gotten
% from raw data);
CorrSong_PBS_pairs=[];
CorrSong_MUSC_pairs=[];

CorrSong_PBS_pairs_OldVer = [];
CorrSong_MUSC_pairs_OldVer = [];
                    CorrMotif_PBS_pairs_OldVer = [];
                    CorrMotif_MUSC_pairs_OldVer = [];

CorrMotif_PBS_pairs=[];
CorrMotif_MUSC_pairs=[];
PairedAcousticDist=[];
PairedAcousticDist_MUSC=[];
IsSameSyl = [];

IsSameMotif = [];
NumSylsInBetween = [];

NumSylsSharedInContext = [];


exptcount=0;
counttmp=1;
if NullControl_splitday==0
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
    
    
    for i=1:NumBirds
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        exptstoplot = 1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        if onlyUseFirstExpt==1
            if max(exptstoplot)>1
                % ================= MODIFY TO TAKE FIRST EXPT WITH MAXIMUM
                % NUMBER OF SYLS
                %             numexperiments=1;
                
                NumExptList = {'pu11wh87', 3, ...
                    'gr41gr90', 1, ...
                    'rd23gr89', 2, ...
                    'rd28pu64', 1, ...
                    'bk34bk68', 1, ...
                    }; % earliest expt with max num syls presnet.
                
                ind = find(strcmp(NumExptList, birdname));
                exptstoplot = NumExptList{ind+1};
                
            end
            
        end
        
        for ii=exptstoplot
            exptcount=exptcount+1;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            % ================================= LOOK AT DEFINED MOTIFS
            % ----------- compare regexp motifs and hand entered
            % motifs
            disp('=================');
            nmHAND = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.motifs);
            nmREG = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList);
            if nmHAND ~= nmREG
                disp(['[NOTE:!!! : num motifs (hand -- regexp) = ' num2str(nmHAND) '-' num2str(nmREG)])
            end
            
            if (0) % === TURN ON to compare regexp and motifs
                for nm = 1:nmHAND
                    disp('-----');
                    disp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.motifs{nm});
                    disp([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList{nm}]);
                end
            end
            
            % =========== NOTE TO SELF: motifs match for regexp and hand (a few rare cases
            % one has more syls than the other, but the orders are all
            % identical. also the number of motifs is identical.
            
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
                
                syl1_single = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
                SingleSylAll = [SingleSylAll syl1_single];
                
                
                % ======================= LEARNING RELATED STUFF
                learnhz = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                
                targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                learnhz_targ = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
                
                LearnHzAll = [LearnHzAll learnhz];
                LearnHzTargAll = [LearnHzTargAll learnhz_targ];

                istarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                IsTargAll = [IsTargAll istarg];
                
                targlearndir = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                TargLearnDir = [TargLearnDir targlearndir];

                distfromtarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
                DistFromTarg = [DistFromTarg distfromtarg];
                
                % ---- extract paired stuff
                for jj=j+1:length(SylsUnique)
                    syl2=SylsUnique{jj};
                    
                    syl2_single = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).single_syl;
                    
                    % ============== if same syl, how many syls back?
                    
                    
                    
                    % --- sanity check - so that matches next section
                    % paired code
                    acoustictmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_zscore_mean;
                    acousticdist=sqrt(sum((acoustic_vec_PBS-acoustictmp).^2));
                    
                    acoustictmp_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl2).fv_baseline_zscore_mean;
                    acousticdist_MUSC=sqrt(sum((acoustic_vec_MUSC-acoustictmp_musc).^2));
                    
                    corrsongPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                    corrsongMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                    
                    
                    CorrSong_PBS_pairs_OldVer = [CorrSong_PBS_pairs_OldVer corrsongPBS];
                    CorrSong_MUSC_pairs_OldVer = [CorrSong_MUSC_pairs_OldVer corrsongMUSC];
                    
                    
                    % ******** MOTIF BY MOTIF CORRELATION
                    corrmotifPBS=nan;
                    corrmotifMUSC=nan;
                    
                    if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS, 'motif_by_motif');
                        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, syl2);
                            
                            corrmotifPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                            corrmotifMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                        end
                    end
                    
                    CorrMotif_PBS_pairs_OldVer = [CorrMotif_PBS_pairs_OldVer corrmotifPBS];
                    CorrMotif_MUSC_pairs_OldVer = [CorrMotif_MUSC_pairs_OldVer corrmotifMUSC];
                    
                    %% sanity check of corr === REXTRACT ORIGINAL DATA USED TO CALCULATE CORR
                    
                    %===================== 1) PBS
                    % ----- confirm that this correaltion value is taken
                    % from this older structure
                    syl1Pos = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.ind_of_this_syl;
                    syl2Pos = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).CORRELATIONS.song_by_song.ind_of_this_syl;
                    
                    if ~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2))
                        assert(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2) == ...
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.RhoMat(syl1Pos, syl2Pos), 'problem?');
                        
                        % ------ extract original raw FF values
                        ffvals1ORIG = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled(:, syl1Pos);
                        ffvals2ORIG = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled(:, syl2Pos);
                        ffvals1_nonanORIG = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled_NanSongsRemoved(:, syl1Pos);
                        ffvals2_nonanORIG = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled_NanSongsRemoved(:, syl2Pos);
                        
                        
                        tvalsORIG = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.Song_datenums;
                        %                     tvals2ORIG = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.Song_datenums;
                        
                        assert(abs(corr(ffvals1_nonanORIG, ffvals2_nonanORIG) - corrsongPBS) < 0.01);
                        
                        % *************** MOTIF BY MOTIF
                        motifnum1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_motifnum;
                        sylpos1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_PosInMotif_thissyl;
                        
                        motifnum2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_motifnum;
                        sylpos2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_PosInMotif_thissyl;
                        
                        if motifnum1~=motifnum2;
                            % then not on same motif, can't calc motif corr
                            doMotifCorr=0;
                        else
                            doMotifCorr=1;
                        end
                            
                        if doMotifCorr==1
                        % =========== PBS
                        ffvalsMOTIF = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{motifnum1}.sub_class{1}.FFvals(:,[sylpos1 sylpos2]);
                        tvalsMOTIF = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{motifnum1}.sub_class{1}.Tvals;
                        
                        % ========= MUSC
                        ffvalsMOTIF_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.data_ParsedIntoSubclasses{motifnum1}.sub_class{1}.FFvals(:,[sylpos1 sylpos2]);
                        tvalsMOTIF_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.data_ParsedIntoSubclasses{motifnum1}.sub_class{1}.Tvals;
                        
                        % --- sanity check
%                         ffvals1MOTIF = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data{1}.data_WithOutlier{1};
                        else
                            
                        end
                        
                        
                        %===================== 1) musc
                        % ----- confirm that this correaltion value is taken
                        % from this older structure
                        syl1Pos = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.ind_of_this_syl;
                        syl2Pos = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl2).CORRELATIONS.song_by_song.ind_of_this_syl;
                        
                        assert(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2) == ...
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.RhoMat(syl1Pos, syl2Pos), 'problem?');
                        
                        % ------ extract original raw FF values
                        ffvals1ORIG_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled(:, syl1Pos);
                        ffvals2ORIG_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled(:, syl2Pos);
                        ffvals1_nonanORIG_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled_NanSongsRemoved(:, syl1Pos);
                        ffvals2_nonanORIG_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled_NanSongsRemoved(:, syl2Pos);
                        
                        tvalsORIG_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.Song_datenums;
                        %                     tvals2ORIG_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.Correlations_Across_Classes.SONG_BY_SONG.Song_datenums(:, syl2Pos);
                        
                        assert(abs(corr(ffvals1_nonanORIG_MUSC, ffvals2_nonanORIG_MUSC) -corrsongMUSC) < 0.01);
                        
                        
                        
                        
                        % ===================== ONLY KEEP DAYS WITH BOTH
                        % PBS AND MUSC
                        if OnlyKeepDaysWithMUSC ==1
                            
                            % =========== SONG BY SONG
                            dayswithmusc = unique(floor(tvalsORIG_MUSC));
                            indstokeepPBS = ismember(floor(tvalsORIG), dayswithmusc);
                            
                            ffvals1ORIG = ffvals1ORIG(indstokeepPBS);
                            ffvals2ORIG = ffvals2ORIG(indstokeepPBS);
                            tvalsORIG = tvalsORIG(indstokeepPBS);
                            
                            % =========== MOTIF BY MOTIF
                            dayswithmusc = unique(floor(tvalsMOTIF_MUSC));
                            indstokeepPBS = ismember(floor(tvalsMOTIF), dayswithmusc);
                            
                            ffvalsMOTIF = ffvalsMOTIF(indstokeepPBS, :);
                            tvalsMOTIF = tvalsMOTIF(indstokeepPBS);
                        end
                        
                        % ################################## SUBTRACT DAY MEANS
                        if CorrSubtractDayMean==1
                            % ================================== PBS
                            tdays = floor(tvalsORIG);
                            ff1 = ffvals1ORIG;
                            ff2 = ffvals2ORIG;
                            
                            % -- get mean across days
                            ff1_daymeans = grpstats(ff1, tdays, {'mean'});
                            ff2_daymeans = grpstats(ff2, tdays, {'mean'});
                            
                            % -- subtract daymean from each value
                            ff1_daymeans_trials = ones(size(ff1));
                            eachday = unique(tdays);
                            for day = eachday
                                ff1_daymeans_trials(tdays ==day) = ff1_daymeans(eachday==day);
                            end
                            ff1_minusday = ff1 - ff1_daymeans_trials;
                            
                            % -- subtract daymean from each value
                            ff2_daymeans_trials = ones(size(ff2));
                            eachday = unique(tdays);
                            for day = eachday
                                ff2_daymeans_trials(tdays ==day) = ff2_daymeans(eachday==day);
                            end
                            ff2_minusday = ff2 - ff2_daymeans_trials;
                            
                            ffvals1ORIG = ff1_minusday;
                            ffvals2ORIG = ff2_minusday;
                            
                            % ================================== MUSC
                            tdays = floor(tvalsORIG_MUSC);
                            ff1 = ffvals1ORIG_MUSC;
                            ff2 = ffvals2ORIG_MUSC;
                            
                            % -- get mean across days
                            ff1_daymeans = grpstats(ff1, tdays, {'mean'});
                            ff2_daymeans = grpstats(ff2, tdays, {'mean'});
                            
                            % -- subtract daymean from each value
                            ff1_daymeans_trials = ones(size(ff1));
                            eachday = unique(tdays);
                            for day = eachday
                                ff1_daymeans_trials(tdays ==day) = ff1_daymeans(eachday==day);
                            end
                            ff1_minusday = ff1 - ff1_daymeans_trials;
                            
                            % -- subtract daymean from each value
                            ff2_daymeans_trials = ones(size(ff2));
                            eachday = unique(tdays);
                            for day = eachday
                                ff2_daymeans_trials(tdays ==day) = ff2_daymeans(eachday==day);
                            end
                            ff2_minusday = ff2 - ff2_daymeans_trials;
                            
                            ffvals1ORIG_MUSC = ff1_minusday;
                            ffvals2ORIG_MUSC = ff2_minusday;
                            
                            % ********************** MOTIF BY MOTIF
                            % ================================== PBS
                            tdays = floor(tvalsMOTIF);
                            ff1 = ffvalsMOTIF;
                            
                            % -- get mean across days
                            ff1_daymeans = grpstats(ff1, tdays, {'mean'});
                            
                            % -- subtract daymean from each value
                            ff1_daymeans_trials = ones(size(ff1));
                            eachday = unique(tdays);
                            for day = eachday'
                                ff1_daymeans_trials(tdays==day,:) = repmat(ff1_daymeans(eachday==day,:), sum(tdays==day), 1);
                            end
                            ffvalsMOTIF = ff1 - ff1_daymeans_trials;
                            
                            % ================================== PBS
                            tdays = floor(tvalsMOTIF_MUSC);
                            ff1 = ffvalsMOTIF_MUSC;
                            
                            % -- get mean across days
                            ff1_daymeans = grpstats(ff1, tdays, {'mean'});
                            
                            % -- subtract daymean from each value
                            ff1_daymeans_trials = ones(size(ff1));
                            eachday = unique(tdays);
                            for day = eachday'
                                ff1_daymeans_trials(tdays==day,:) = repmat(ff1_daymeans(eachday==day,:), sum(tdays==day), 1);
                            end
                            ffvalsMOTIF_MUSC = ff1 - ff1_daymeans_trials;

                        end
                        
                        % ================== RECALC
                        if reCalcOldMethodKeepingNanSongs ==1
                            % =========== then recalculate ...
                            % ======================= SONG BY SONG
                            indtmp = ~isnan(ffvals1ORIG) & ~isnan(ffvals2ORIG);
                            corrsongPBS = corr(ffvals1ORIG(indtmp), ffvals2ORIG(indtmp));
                            
                            indtmp = ~isnan(ffvals1ORIG_MUSC) & ~isnan(ffvals2ORIG_MUSC);
                            corrsongMUSC = corr(ffvals1ORIG_MUSC(indtmp), ffvals2ORIG_MUSC(indtmp));
                            
                            % ======================== MOTIF
                            % ------ PBS
                            indtmp = ~isnan(ffvalsMOTIF(:,1)) & ~isnan(ffvalsMOTIF(:,2));
                            corrmotifPBS = corr(ffvalsMOTIF(indtmp,1), ffvalsMOTIF(indtmp,2));
                            
                            % ------ MUSC
                            indtmp = ~isnan(ffvalsMOTIF_MUSC(:,1)) & ~isnan(ffvalsMOTIF_MUSC(:,2));
                            corrmotifMUSC = corr(ffvalsMOTIF_MUSC(indtmp,1), ffvalsMOTIF_MUSC(indtmp,2));
                            
                            
                        end
                        
                        
                    end
                    %% === get day to day FF
                    daysMUSC = find(~cellfun('isempty', ...
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds));
                    dayFirstWN = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                    daysMUSC = daysMUSC(daysMUSC<dayFirstWN);
                    
                    daysBase = 1:dayFirstWN-1;
                    
                    syldat = struct;
                    listofsyls = {syl, syl2};
                    for sylthis = listofsyls;
                        tvals_PBS = [];
                        ffvals_PBS = [];
                        tvals_MUSC = [];
                        ffvals_MUSC = [];
                        
                        ffdiff = []; % one value for each day
                        
                        for ddd = daysMUSC
                            % ============= MUSC
                            ff2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis{1}).FFvals_WithinTimeWindow{ddd};
                            tt = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis{1}).Tvals_WithinTimeWindow{ddd};
                            
                            if isempty(ff2)
                                continue
                            end
                            ffvals_MUSC = [ffvals_MUSC ff2];
                            tvals_MUSC = [tvals_MUSC tt];
                            
                            % =========== PBS
                            ff1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis{1}).FFvals_WithinTimeWindow{ddd};
                            tt = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis{1}).Tvals_WithinTimeWindow{ddd};
                            
                            ffvals_PBS = [ffvals_PBS ff1];
                            tvals_PBS = [tvals_PBS tt];
                            
                            
                            % --- difference
                            ffdiff = [ffdiff mean(ff2)-mean(ff1)];
                        end
                        
                        
                        
                        if (1)
                            % -------------- use all days for PBS (not
                            % just inactivation days)
                            % -- NOTE: THIS ONLY FOR CORR (diff will
                            % still by mean of day means
                            tvals_PBS = [];
                            ffvals_PBS = [];
                            daysBase = 1:dayFirstWN-1;
                            
                            for ddd = daysBase
                                
                                % =========== PBS
                                ff1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis{1}).FFvals_WithinTimeWindow{ddd};
                                tt = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis{1}).Tvals_WithinTimeWindow{ddd};
                                
                                ffvals_PBS = [ffvals_PBS ff1];
                                tvals_PBS = [tvals_PBS tt];
                            end
                        end
                        
                        
                        % ---- save to syl dat
                        syldat.(sylthis{1}).tvals_PBS = tvals_PBS;
                        syldat.(sylthis{1}).ffvals_PBS = ffvals_PBS;
                        syldat.(sylthis{1}).tvals_MUSC = tvals_MUSC;
                        syldat.(sylthis{1}).ffvals_MUSC = ffvals_MUSC;
                        syldat.(sylthis{1}).ffmean_MUSCminusPBS = ffdiff;
                        
                        
                        % ########################### get all days
                        tvals_PBS_allBase = [];
                        ffvals_PBS_allBase = [];
                        tvals_MUSC_allBase = [];
                        ffvals_MUSC_allBase = [];
                        
                        for ddd = daysBase
                            % ============= MUSC
                            ff2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis{1}).FFvals_WithinTimeWindow{ddd};
                            tt = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis{1}).Tvals_WithinTimeWindow{ddd};
                            
                            %                             if isempty(ff2)
                            %                                 continue
                            %                             end
                            ffvals_MUSC_allBase = [ffvals_MUSC_allBase ff2];
                            tvals_MUSC_allBase = [tvals_MUSC_allBase tt];
                            
                            % =========== PBS
                            ff1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis{1}).FFvals_WithinTimeWindow{ddd};
                            tt = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis{1}).Tvals_WithinTimeWindow{ddd};
                            
                            ffvals_PBS_allBase = [ffvals_PBS_allBase ff1];
                            tvals_PBS_allBase = [tvals_PBS_allBase tt];
                            
                        end
                        
                        syldat.(sylthis{1}).tvals_PBS_allBase = tvals_PBS_allBase;
                        syldat.(sylthis{1}).ffvals_PBS_allBase = ffvals_PBS_allBase;
                        syldat.(sylthis{1}).tvals_MUSC_allBase = tvals_MUSC_allBase;
                        syldat.(sylthis{1}).ffvals_MUSC_allBase = ffvals_MUSC_allBase;
                    end
                    
                    % ############################## 1) RECALCULATE
                    % CHANGE IN FF
                    % --------------- ONLY FOR SYL 1
                    baseline_musc_minus_pbs = mean(syldat.(syl).ffmean_MUSCminusPBS);
                    
                    % ############################# SAVE DAY BY DAY (JUST FOR SYL 1)
                    if jj == j+1
                        % -- [this gets 1:N-1 syls] then save first syl, do not save for jj
                        % greater than j+1 since that will be redundant
                        PitchTvalsPBS = [PitchTvalsPBS syldat.(syl).tvals_PBS_allBase];
                        PitchFFvalsPBS = [PitchFFvalsPBS syldat.(syl).ffvals_PBS_allBase];
                        
                        PitchTvalsMUSC = [PitchTvalsMUSC syldat.(syl).tvals_MUSC_allBase];
                        PitchFFvalsMUSC = [PitchFFvalsMUSC syldat.(syl).ffvals_MUSC_allBase];
                    end
                    if j==length(SylsUnique)-1 & jj==j+1
                        % -- [this gets syl number N]
                        PitchTvalsPBS = [PitchTvalsPBS syldat.(syl2).tvals_PBS_allBase];
                        PitchFFvalsPBS = [PitchFFvalsPBS syldat.(syl2).ffvals_PBS_allBase];
                        
                        PitchTvalsMUSC = [PitchTvalsMUSC syldat.(syl2).tvals_MUSC_allBase];
                        PitchFFvalsMUSC = [PitchFFvalsMUSC syldat.(syl2).ffvals_MUSC_allBase];
                    end
                    
                    
                    
                    if recalcValues ==1
                        %% RECALCULATE CORRELATIONS,
                        % ALSO GETTING CROSS CORRELATION (I.E. AT DIFFERENT
                        % LAGS)
                        
                        % ================= 1) extract baseline time and
                        % ff
                        % --- only collect from days with LMAN inactivation
                        
                        
                        % ############################# 2) RECALC CORR
                        % =================== NEW VERSION
                        
                        % ############################# get song by song ff
                        % 1) --------- syl1, PBS
                        sylthis = syl;
                        tfield = 'tvals_PBS';
                        ffield = 'ffvals_PBS';
                        
                        syldat = fn_ffsongbysong(syldat, sylthis, tfield, ffield);
                        
                        % 2) --------- syl1, MUSC
                        sylthis = syl;
                        tfield = 'tvals_MUSC';
                        ffield = 'ffvals_MUSC';
                        
                        syldat = fn_ffsongbysong(syldat, sylthis, tfield, ffield);
                        
                        % 3) --------- syl2, PBS
                        sylthis = syl2;
                        tfield = 'tvals_PBS';
                        ffield = 'ffvals_PBS';
                        
                        syldat = fn_ffsongbysong(syldat, sylthis, tfield, ffield);
                        
                        % 4) --------- syl2, MUSC
                        sylthis = syl2;
                        tfield = 'tvals_MUSC';
                        ffield = 'ffvals_MUSC';
                        
                        syldat = fn_ffsongbysong(syldat, sylthis, tfield, ffield);
                        
                        %% sanity check - compare reextracted to original data
                        if rand<0
                            hsplots = [];
                            lt_figure; hold on;
                            % ################# OLD METHOD
                            hsplot = lt_subplot(2,1,1); hold on;
                            hsplots = [hsplots hsplot];
                            title('old method');
                            % ========== PBS
                            plot(tvalsORIG, ffvals1ORIG, 'ok');
                            plot(tvalsORIG, ffvals2ORIG, 'sk');
                            
                            % =========== MUSC
                            plot(tvalsORIG_MUSC, ffvals1ORIG_MUSC, 'or');
                            plot(tvalsORIG_MUSC, ffvals2ORIG_MUSC, 'sr');
                            
                            lt_plot_annotation(1, ['rho=' num2str(corrsongPBS)]);
                            lt_plot_annotation(2, ['rho=' num2str(corrsongMUSC)], 'r');
                            
                            
                            % ################## NEW METHOD
                            hsplot = lt_subplot(2,1,2); hold on;
                            hsplots = [hsplots hsplot];
                            title('new method');
                            % ========== PBS
                            plot(syldat.(syl).tvals_PBS_song, syldat.(syl).ffvals_PBS_song, 'ok');
                            plot(syldat.(syl2).tvals_PBS_song, syldat.(syl2).ffvals_PBS_song, 'sk');
                            
                            % =========== MUSC
                            plot(syldat.(syl).tvals_MUSC_song, syldat.(syl).ffvals_MUSC_song, 'or');
                            plot(syldat.(syl2).tvals_MUSC_song, syldat.(syl2).ffvals_MUSC_song, 'sr');
                            
                            % ==== annotate correlation
                            % 1) --------------- PBS
                            tfield = 'tvals_PBS_song';
                            ffield = 'ffvals_PBS_song';
                            corrsongPBSTMP = fn_getsongcorr(syldat, syl, syl2, tfield, ffield);
                            
                            % 2) ----------------- MUSC
                            tfield = 'tvals_MUSC_song';
                            ffield = 'ffvals_MUSC_song';
                            corrsongMUSCTMP = fn_getsongcorr(syldat, syl, syl2, tfield, ffield);
                            lt_plot_annotation(1, ['rho=' num2str(corrsongPBSTMP)]);
                            lt_plot_annotation(2, ['rho=' num2str(corrsongMUSCTMP)], 'r');
                            
                            
                            
                            linkaxes(hsplots, 'xy');
                            
                            pause; close all;
                            
                        end
                        
                        
                        %% --- do day mean subtraction
                        
                        % ####################  subtract day means
                        % 1) --------- syl1, PBS
                        sylthis = syl;
                        tfield = 'tvals_PBS_song';
                        ffield = 'ffvals_PBS_song';
                        
                        syldat = fn_subtractDayMean(syldat, sylthis, tfield, ffield);
                        
                        % 2) --------- syl1, MUSC
                        sylthis = syl;
                        tfield = 'tvals_MUSC_song';
                        ffield = 'ffvals_MUSC_song';
                        
                        syldat = fn_subtractDayMean(syldat, sylthis, tfield, ffield);
                        
                        % 3) --------- syl2, PBS
                        sylthis = syl2;
                        tfield = 'tvals_PBS_song';
                        ffield = 'ffvals_PBS_song';
                        
                        syldat = fn_subtractDayMean(syldat, sylthis, tfield, ffield);
                        
                        % 4) --------- syl2, MUSC
                        sylthis = syl2;
                        tfield = 'tvals_MUSC_song';
                        ffield = 'ffvals_MUSC_song';
                        
                        syldat = fn_subtractDayMean(syldat, sylthis, tfield, ffield);
                        
                        %%  ############ CALCULATE CORR
                        
                        corrsongPBS = nan;
                        corrsongMUSC =nan;
                        corrmotifPBS=nan;
                        corrmotifMUSC=nan;
                        
                        if CorrSubtractDayMean ==1
                            
                            
                            % 1) --------------- PBS
                            tfield = 'tvals_PBS_song';
                            ffield = 'ffvals_PBS_song_minusday';
                            corrsongPBS = fn_getsongcorr(syldat, syl, syl2, tfield, ffield);
                            
                            % 2) ----------------- MUSC
                            tfield = 'tvals_MUSC_song';
                            ffield = 'ffvals_MUSC_song_minusday';
                            corrsongMUSC = fn_getsongcorr(syldat, syl, syl2, tfield, ffield);
                        else
                            % ----- don't do day subtraction
                            
                            % 1) --------------- PBS
                            tfield = 'tvals_PBS_song';
                            ffield = 'ffvals_PBS_song';
                            corrsongPBS = fn_getsongcorr(syldat, syl, syl2, tfield, ffield);
                            
                            % 2) ----------------- MUSC
                            tfield = 'tvals_MUSC_song';
                            ffield = 'ffvals_MUSC_song';
                            corrsongMUSC = fn_getsongcorr(syldat, syl, syl2, tfield, ffield);
                            
                        end
                        
                        
                    end
                    
                    % === are they same type pair? (based on single syl)
                    samepair = strcmp(syl1_single, syl2_single);
                    
                    SylID1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl);
                    SylID2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2);
                    %                     motifs1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.
                    
                    % ======= if same syl, how many syls shared in
                    % context?
                    sylsshared = nan;
                        if strcmp(SylID1.preceding_syl, SylID2.preceding_syl)==0
                            % -- no syls shared
                            sylsshared=0;
                        elseif strcmp(SylID1.preceding_syl, SylID2.preceding_syl)==1 ...
                                & strcmp(SylID1.two_syl_back, SylID2.two_syl_back)==0
                            % then only one syl shared
                            sylsshared=1;
                        elseif strcmp(SylID1.preceding_syl, SylID2.preceding_syl)==1 ...
                                & strcmp(SylID1.two_syl_back, SylID2.two_syl_back)==1
                            % then two syls shared
                            sylsshared=2;
                        elseif strcmp(SylID1.preceding_syl, SylID2.preceding_syl)==1 ...
                                & strcmp(SylID1.two_syl_back, SylID2.two_syl_back)==1 ...
                                & strcmp(SylID1.three_syl_back, SylID2.three_syl_back)==1
                            sylsshared=3;
                            
                            % === figure out if actually 3 syls shared
                            %                            SylID1
                        end
                    
                    if doMotifCorr==0
                        corrmotifPBS = nan;
                        corrmotifMUSC = nan;
                    end
                       
                    % ======= are they on same motif? if so how far apart?
                    issamemotif = nan;
                    numsylsbetween = nan;
                    if SylID1.motif_num == SylID2.motif_num;
                        issamemotif=1;
                        numsylsbetween = abs(SylID1.regexp_PosInMotif_thissyl - SylID2.regexp_PosInMotif_thissyl) - 1;
                    else
                        issamemotif=0;
                        numsylsbetween = nan;
                    end
                    
                    if isempty(numsylsbetween)
                        numsylsbetween = nan;
                    end
                    
                    % ===== OUTPUT
                    if UseMotifByMotifCorr ==1
                        corrsongPBS = corrmotifPBS;
                        corrsongMUSC = corrmotifMUSC;
                    end
                    CorrMotif_PBS_pairs=[CorrMotif_PBS_pairs corrmotifPBS];
                    CorrMotif_MUSC_pairs=[CorrMotif_MUSC_pairs corrmotifMUSC];
                    CorrSong_PBS_pairs=[CorrSong_PBS_pairs corrsongPBS];
                    CorrSong_MUSC_pairs=[CorrSong_MUSC_pairs corrsongMUSC];
                    PairedAcousticDist=[PairedAcousticDist acousticdist];
                    counttmp=[counttmp counttmp(end)+1];
                    
                    PairedAcousticDist_MUSC=[PairedAcousticDist_MUSC acousticdist_MUSC];
                    IsSameSyl = [IsSameSyl samepair];
                    
                    IsSameMotif = [IsSameMotif issamemotif];
                    NumSylsInBetween = [NumSylsInBetween numsylsbetween];
                    
                    NumSylsSharedInContext = [NumSylsSharedInContext sylsshared];
                end
                
                
                % ======= COLLECT [single syl]
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
    warned = 0;
    
    if UseLMANDatForSplit==1
        filter = 'LMAN';
        [SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
    end
    
    NumBirds=length(SeqDepPitch_AcrossBirds.birds);
    % --- then get mock data
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        %         numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        exptstoplot = 1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        if onlyUseFirstExpt==1 & UseLMANDatForSplit==1
            if max(exptstoplot)>1
                % ================= MODIFY TO TAKE FIRST EXPT WITH MAXIMUM
                % NUMBER OF SYLS
                %             numexperiments=1;
                
                NumExptList = {'pu11wh87', 3, ...
                    'gr41gr90', 1, ...
                    'rd23gr89', 2, ...
                    'rd28pu64', 1, ...
                    'bk34bk68', 1, ...
                    }; % earliest expt with max num syls presnet.
                
                ind = find(strcmp(NumExptList, birdname));
                exptstoplot = NumExptList{ind+1};
                
            end
            
        elseif onlyUseFirstExpt==1 & UseLMANDatForSplit==0 & warned==0
            warned=1;
            disp('IMPROTANT: have not yet implemented onlyUseFirstExpt for split analyses');
            pause;
        end
        
        %         for ii=1:numexperiments;
        for ii=exptstoplot
            
            % === skip if is LMAN experiment
            if UseLMANDatForSplit==1
            else
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                    continue
                end
            end
            
            
            exptcount=exptcount+1;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            for j=1:length(SylsUnique)
                syl=SylsUnique{j};
                
                % determine which tvals are mock "PBS" and which are mock
                % "MUSC"
                if UseLMANDatForSplit==0
                    baseline_tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals;
                    baseline_ffvals_all=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF;
                else
                    % --- only collect from days without LMAN inactivation
                    daysPBS = find(cellfun('isempty', ...
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds));
                    dayFirstWN = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                    daysPBS = daysPBS(daysPBS<dayFirstWN);
                    
                    baseline_tvals = [];
                    baseline_ffvals_all = [];
                    for ddd = daysPBS
                        fvals = cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{ddd});
                        tvals = cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{ddd});
                        
                        baseline_tvals = [baseline_tvals tvals];
                        baseline_ffvals_all = [baseline_ffvals_all fvals];
                    end
                end
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
                baseline_pbs=mean(baseline_ffvals_all(inds_PBS));
                pitchstd_pbs=std(baseline_ffvals_all(inds_PBS));
                
                baseline_musc=mean(baseline_ffvals_all(inds_MUSC));
                
                baseline_musc_minus_pbs=baseline_musc-baseline_pbs;
                
                acoustic_vec_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore_mean;
                
                
                singleSyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
                SingleSylAll = [SingleSylAll singleSyl];
                
                
                
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
IsSameSyl = logical(IsSameSyl);


%% ================== SAVE INTO OUTPUT STRUCT

OUTSTRUCT.singlesyls.AcousVec_PBS=AcousVec_PBS;
OUTSTRUCT.singlesyls.Pitch_MUSC_minus_PBS=Pitch_MUSC_minus_PBS;
OUTSTRUCT.singlesyls.Pitch_MUSC_minus_PBS_zscore = Pitch_MUSC_minus_PBS_zscore;
OUTSTRUCT.singlesyls.Pitch_PBSraw=Pitch_PBSraw;
OUTSTRUCT.singlesyls.PitchSTD_PBS=PitchSTD_PBS;
OUTSTRUCT.singlesyls.ExptCounter=ExptCounter;
OUTSTRUCT.singlesyls.Birdnames=Birdnames;
OUTSTRUCT.singlesyls.Exptnames=Exptnames;
OUTSTRUCT.singlesyls.Syllables=Syllables;
OUTSTRUCT.singlesyls.SingleSylAll = SingleSylAll;
OUTSTRUCT.singlesyls.PitchTvalsPBS = PitchTvalsPBS;
OUTSTRUCT.singlesyls.PitchTvalsMUSC = PitchTvalsMUSC;
OUTSTRUCT.singlesyls.PitchFFvalsPBS = PitchFFvalsPBS;
OUTSTRUCT.singlesyls.PitchFFvalsMUSC = PitchFFvalsMUSC;
                OUTSTRUCT.singlesyls.LearnHzAll = LearnHzAll;
                OUTSTRUCT.singlesyls.LearnHzTargAll = LearnHzTargAll;
                OUTSTRUCT.singlesyls.IsTargAll = IsTargAll;
OUTSTRUCT.singlesyls.TargLearnDir = TargLearnDir;
             OUTSTRUCT.singlesyls.DistFromTarg = DistFromTarg;
             
             
OUTSTRUCT.pairedsyls.CorrSong_PBS_pairs=CorrSong_PBS_pairs;
OUTSTRUCT.pairedsyls.CorrSong_MUSC_pairs=CorrSong_MUSC_pairs;
OUTSTRUCT.pairedsyls.CorrSong_PBS_pairs_OldVer =CorrSong_PBS_pairs_OldVer;
OUTSTRUCT.pairedsyls.CorrSong_MUSC_pairs_OldVer = CorrSong_MUSC_pairs_OldVer;
OUTSTRUCT.pairedsyls.CorrMotif_PBS_pairs=CorrMotif_PBS_pairs;
OUTSTRUCT.pairedsyls.CorrMotif_MUSC_pairs=CorrMotif_MUSC_pairs;
OUTSTRUCT.pairedsyls.PairedAcousticDist=PairedAcousticDist;
OUTSTRUCT.pairedsyls.PairedAcousticDist_MUSC=PairedAcousticDist_MUSC;
OUTSTRUCT.pairedsyls.IsSameSyl = IsSameSyl;
OUTSTRUCT.pairedsyls.IsSameMotif = IsSameMotif;
OUTSTRUCT.pairedsyls.NumSylsInBetween = NumSylsInBetween;
OUTSTRUCT.pairedsyls.NumSylsSharedInContext = NumSylsSharedInContext;


%% ###### compare new vs. old method correlations
    if reCalcOldMethodKeepingNanSongs ==1
        lt_figure; hold on
    
    % ==== pbs
    lt_subplot(4,1,1); hold on;
    title('pbs (song)');
    xlabel('using regexp')
    ylabel('using new ver')
    plot(CorrSong_PBS_pairs_OldVer, CorrSong_PBS_pairs, 'ok');
    
    
    % ==== musc
    lt_subplot(4,1,2); hold on;
    title('musc (song)');
    plot(CorrSong_MUSC_pairs_OldVer, CorrSong_MUSC_pairs, 'ok');
    
    
    lt_subplot(4,1,3); hold on;
    title('PBS (motif)');
    plot(CorrMotif_PBS_pairs_OldVer, CorrMotif_PBS_pairs, 'ok');

    lt_subplot(4,1,4); hold on;
    title('musc (motif)');
    plot(CorrMotif_MUSC_pairs_OldVer, CorrMotif_MUSC_pairs, 'ok');
    
    end

%% =============== FIRST, collect data (for all pairs)
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
Pairs_OriginalInds = [];

% for split analysis only
IsSameSylv2 = [];

% --- for each experiment, get all pairwise 1) acoustic distance, and 2)
% similarity of shift
NumExpts=max(ExptCounter);

for i=1:NumExpts
    Inds=find(ExptCounter==i);
    
    % ---- get all pairwise stuff
    for j=1:length(Inds)
        ind1=Inds(j);
        
        for jj=j+1:length(Inds)
            ind2=Inds(jj);
            
            acousticDist=sqrt(sum((AcousVec_PBS(:, ind1)-AcousVec_PBS(:, ind2)).^2));
            muscEffectDiff=abs(Pitch_MUSC_minus_PBS(ind1)-Pitch_MUSC_minus_PBS(ind2));
            similarDirOfShift=sign(Pitch_MUSC_minus_PBS(ind1))==sign(Pitch_MUSC_minus_PBS(ind2));
            
            separationPBS=Pitch_PBSraw(ind1)-Pitch_PBSraw(ind2);
            separationMUSC=(Pitch_PBSraw(ind1)+Pitch_MUSC_minus_PBS(ind1)) - (Pitch_PBSraw(ind2)+Pitch_MUSC_minus_PBS(ind2));
            changeInAbsSeparationDuringMUSC=abs(separationMUSC)-abs(separationPBS);
            
            IsSameSylv2 = [IsSameSylv2 strcmp(SingleSylAll{ind1}, SingleSylAll{ind2})];
            
            % ==== collect pairwise stats
            PitchShiftBothSyls_paired=[PitchShiftBothSyls_paired [Pitch_MUSC_minus_PBS(ind1) Pitch_MUSC_minus_PBS(ind2)]];
            AcousticDist_pairs=[AcousticDist_pairs acousticDist];
            DiffInPitchShiftCausedByMusc_pairs=[DiffInPitchShiftCausedByMusc_pairs muscEffectDiff];
            similarDirOfShift_pairs=[similarDirOfShift_pairs similarDirOfShift];
            changeInAbsSeparationDuringMUSC_pairs=[changeInAbsSeparationDuringMUSC_pairs changeInAbsSeparationDuringMUSC];
            separationPBS_paired=[separationPBS_paired separationPBS];
            separationMUSC_paired=[separationMUSC_paired separationMUSC];
            Pairs_OriginalInds = [Pairs_OriginalInds; [ind1 ind2]];
            
            %            disp([num2str(ind1) ' - '  num2str(ind2)])
            
        end
    end
end

if NullControl_splitday==0
    assert(all(IsSameSyl == IsSameSylv2), 'problem - not same...');
else
    IsSameSyl = logical(IsSameSylv2);
end

% == confirm that things match from previuos code
if NullControl_splitday==0
    lt_figure; hold on;
    title('checking code match - must be 45deg [if not, then cant trust paired measures]');
    lt_plot_45degScatter(AcousticDist_pairs, PairedAcousticDist);
end

%% =========================== OUTPUT DATA

OUTSTRUCT.pairedsyls.AcousticDist_pairs=AcousticDist_pairs;
OUTSTRUCT.pairedsyls.DiffInPitchShiftCausedByMusc_pairs=DiffInPitchShiftCausedByMusc_pairs;
OUTSTRUCT.pairedsyls.similarDirOfShift_pairs=similarDirOfShift_pairs;
OUTSTRUCT.pairedsyls.changeInAbsSeparationDuringMUSC_pairs=changeInAbsSeparationDuringMUSC_pairs;
OUTSTRUCT.pairedsyls.separationPBS_paired=separationPBS_paired;
OUTSTRUCT.pairedsyls.separationMUSC_paired=separationMUSC_paired;
OUTSTRUCT.pairedsyls.PitchShiftBothSyls_paired=PitchShiftBothSyls_paired;
OUTSTRUCT.pairedsyls.Pairs_OriginalInds = Pairs_OriginalInds;

% for split analysis only
OUTSTRUCT.pairedsyls.IsSameSylv2 = IsSameSylv2;




end

function syldat = fn_ffsongbysong(syldat, sylthis, tfield, ffield)

tvals = syldat.(sylthis).(tfield);
fvals = syldat.(sylthis).(ffield);

fvals_song = grpstats(fvals, tvals, {'mean'});
fvals_song = fvals_song';
tvals_song = unique(tvals);

syldat.(sylthis).([tfield '_song']) = tvals_song;
syldat.(sylthis).([ffield '_song']) = fvals_song;


% -- sanity
if (0)
    lt_figure; hold on;
    plot(tvals, fvals, 'ok');
    lt_plot(tvals_song, fvals_song, {'Color', 'r'});
end



end


function corrsong = fn_getsongcorr(syldat, syl, syl2, tfield, ffield)
tvals1 = syldat.(syl).(tfield);
ffvals1 = syldat.(syl).(ffield);

tvals2 = syldat.(syl2).(tfield);
ffvals2 = syldat.(syl2).(ffield);

% make sure that the tvals are algined (if not,
% then pick out intersect)
[~, inda, indb] = intersect(tvals1, tvals2);
ffvals1 = ffvals1(inda);
ffvals2 = ffvals2(indb);
assert(length(inda)>0.75*length(tvals1), 'why not much intersect?')

corrsong = corr(ffvals1', ffvals2');
assert(~any(isnan(corrsong(:))));
end


function syldat = fn_subtractDayMean(syldat, sylthis, tfield, ffield)

% ---- extract
tdays = floor(syldat.(sylthis).(tfield));
ff = syldat.(sylthis).(ffield);
ffmeans = grpstats(ff, tdays, {'mean'});

% --- get dat sized vector containing day means
ffmeansTrials = ones(size(ff));
eachday = unique(tdays);
for day = eachday
    ffmeansTrials(tdays==day) = ffmeans(eachday==day);
end

% ----- do subtraction
ff_minusday = ff - ffmeansTrials;
syldat.(sylthis).([ffield '_minusday']) = ff_minusday;

% -- sanity
if (0)
    lt_figure; hold on;
    plot(syldat.(sylthis).(tfield), ff, 'ob');
    plot(syldat.(sylthis).(tfield), ffmeansTrials, 'sk');
    pause; close all;
end
end


