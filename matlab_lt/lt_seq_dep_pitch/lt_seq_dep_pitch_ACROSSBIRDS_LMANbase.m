function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbase(SeqDepPitch_AcrossBirds, ...
    PARAMS, same_type_thr, NullControl_splitday, useHandLabForSametype, ...
    onlyUseFirstExpt, recalcValues)
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

% === pairs that cannot be obtained in next section( i.e. must be gotten
% from raw data);
CorrSong_PBS_pairs=[];
CorrSong_MUSC_pairs=[];

CorrSong_PBS_pairs_OldVer = [];
CorrSong_MUSC_pairs_OldVer = [];

CorrMotif_PBS_pairs=[];
CorrMotif_MUSC_pairs=[];
PairedAcousticDist=[];
PairedAcousticDist_MUSC=[];
IsSameSyl = [];

IsSameMotif = [];
NumSylsInBetween = [];

SingleSylAll = {};


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
                
                % ---- extract paired stuff
                for jj=j+1:length(SylsUnique)
                    syl2=SylsUnique{jj};
                    
                    syl2_single = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).single_syl;
                    
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
                    
                    
                    CorrSong_PBS_pairs_OldVer = [CorrSong_PBS_pairs_OldVer corrsongPBS];
                    CorrSong_MUSC_pairs_OldVer = [CorrSong_MUSC_pairs_OldVer corrsongMUSC];
                    
                    
                    
                    if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS, 'motif_by_motif');
                        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, syl2);
                            
                            corrmotifPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                            corrmotifMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                        end
                    end
                    if isempty(corrmotifPBS)
                        keyboard
                    end
                    
                    if recalcValues ==1
                        %% RECALCULATE CORRELATIONS,
                        % ALSO GETTING CROSS CORRELATION (I.E. AT DIFFERENT
                        % LAGS)
                        
                        % ================= 1) extract baseline time and
                        % ff
                        % --- only collect from days with LMAN inactivation
                        daysMUSC = find(~cellfun('isempty', ...
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds));
                        dayFirstWN = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                        daysMUSC = daysMUSC(daysMUSC<dayFirstWN);
                        
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
                            
                            if (0)
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
                        end
                        % ############################## 1) RECALCULATE
                        % CHANGE IN FF
                        % --------------- ONLY FOR SYL 1
                        baseline_musc_minus_pbs = mean(syldat.(syl).ffmean_MUSCminusPBS);
                        
                        
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
                            
                            
                            if (0)
                                % --- do day mean subtraction
                                
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
                                
                                
                                % ############################### CALCULATE
                                % CORR
                                corrsongPBS = nan;
                                corrsongMUSC =nan;
                                corrmotifPBS=nan;
                                corrmotifMUSC=nan;
                                
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
                                
                                % ############################### CALCULATE
                                % CORR
                                corrsongPBS = nan;
                                corrsongMUSC =nan;
                                corrmotifPBS=nan;
                                corrmotifMUSC=nan;
                                
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
                    
                    
                    
                    % ======= are they on same motif? if so how far apart?
                    SylID1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl);
                    SylID2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2);
                    
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


%% ###### compare new vs. old method correlations
if recalcValues==1
   lt_figure; hold on
   
    % ==== pbs
    lt_subplot(2,1,1); hold on;
    title('pbs');
    xlabel('using regexp')
    ylabel('using new ver')
    plot(CorrSong_PBS_pairs_OldVer, CorrSong_PBS_pairs, 'ok');
  
    
    % ==== musc
    lt_subplot(2,1,2); hold on;
    title('musc');
    plot(CorrSong_MUSC_pairs_OldVer, CorrSong_MUSC_pairs, 'ok');
    
    
    
end

%% +++++++++++++ PLOTS

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];



%% ====== DISTRIBUTION OF PBS AND MUSC PITCH
%  ########### 1) distrubitions
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('distr of baseline PBS pitch');

% --- PBS
[~, Xcenters] = lt_plot_histogram(Pitch_PBSraw, '', 1, 1, 0, 1, 'k'); % PBS

% --- MUSC
lt_plot_histogram(Pitch_PBSraw+Pitch_MUSC_minus_PBS, Xcenters, 1, 1, 0, 1, 'r'); % MUSC

% ############## 2) paired
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all syls');

x = [1 2];
y = [Pitch_PBSraw', Pitch_PBSraw'+Pitch_MUSC_minus_PBS'];
plot(x, y, '-ok');
plot(2, y(:,2), 'or');
xlim([0 3]);
xlabel('PBS -- MUSC');


%% =========== DISTRIBUTIONS, SEPARATE FOR EACH EXPERIMENT
figcount=1;
subplotrows=4;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];


numexpts = max(ExptCounter);
for i=1:numexpts
    
    % ############## 2) paired
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    inds = ExptCounter==i;
    
    bname = unique({Birdnames{inds}});
    ename = unique({Exptnames{inds}});
    assert(length([bname ename])==2, 'asfds');
    title([bname '-' ename]);
    
    x = [1 2];
    y = [Pitch_PBSraw(inds)', Pitch_PBSraw(inds)'+Pitch_MUSC_minus_PBS(inds)'];
    sylnames = {Syllables{inds}};
    
    plot(x, y, '-ok');
    plot(2, y(:,2), 'or');
    
    for j=1:length(y)
        lt_plot_text(2.1, y(j,2), sylnames{j}, 'r');
    end
    xlim([0 3]);
    xlabel('PBS -- MUSC');
    
end

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
lt_regress(Pitch_MUSC_minus_PBS, PitchSTD_PBS, 1, 0, 1, 1, 'r', 0);

% --- does shift correlate with baseline std?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('abs(musc shifts) vs. baseline std (across rends and days)');

plot(PitchSTD_PBS, abs(Pitch_MUSC_minus_PBS), 'ok');
lt_regress(abs(Pitch_MUSC_minus_PBS), PitchSTD_PBS, 1, 0, 1, 1, 'r', 0);





%% ################ 2) if you are more similar structurally, do you have more similar musc shift?
%% ###############################################################

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

%% ====== PLOT (scatter of pitch shifts due to musc)


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
% -- randomly switch positions of X and Y (since order doesn't matter)
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
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
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
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
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;


% ################################################################################
% ######################################## same, but using hand label same type
% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (handlab)');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(IsSameSyl);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
color='r';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(IsSameSyl==0);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;


%% ################### [SYL1 = starting higher pitch]

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ######################################## same, but using hand label same type
% === plot all pairwise effect (45 degree) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')

inds=find(IsSameSyl);
XY = cell2mat(PitchShiftBothSyls_paired(inds)');

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));
% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));

% -- annotate cases where MUSC brought separation closer
indsMuscCloser = abs(separationMUSC_paired(inds)) < abs(separationPBS_paired(inds));

% --- PLOT
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
plot(XY(indsMuscCloser, 1), XY(indsMuscCloser,2), 'mo');
lt_plot_zeroline;lt_plot_zeroline_vert;
line([-50 50], [-50 50], 'Color','k');

% --- also plot how much greater pitch of syl 2 is (rel syl 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl1 hz - syl2 hz (base)')

Ysep= abs(ffsyl2-ffsyl1);
plot(XY(:,1), Ysep, 'ok');

% ----------------------------
linkaxes(hsplots, 'x');

% ######################################## diff
% === plot all pairwise effect (45 degree) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('diff [syl1: higher hz]');
color='r';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')

inds=find(IsSameSyl==0);
XY = cell2mat(PitchShiftBothSyls_paired(inds)');

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));
% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));

% -- annotate cases where MUSC brought separation closer
indsMuscCloser = abs(separationMUSC_paired(inds)) < abs(separationPBS_paired(inds));

% --- PLOT
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
plot(XY(indsMuscCloser, 1), XY(indsMuscCloser,2), 'mo');
lt_plot_zeroline;lt_plot_zeroline_vert;
line([-50 50], [-50 50], 'Color','k');
xlim([-55 55]); ylim([-55 55]);

% --- also plot how much greater pitch of syl 2 is (rel syl 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('diff [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl1 hz - syl2 hz (base)')

Ysep= abs(ffsyl2-ffsyl1);
plot(XY(:,1), Ysep, 'ok');

% ----------------------------
linkaxes(hsplots, 'x');


%% #################### DISTRIBUTION OF TYPES OF PAIRED EFFECTS
CategoriesStruct = struct;

% ====== SAME SYL
inds=find(IsSameSyl);

CategoriesStruct.same = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_paired, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 1);

% ====== DIFF SYL
inds=find(IsSameSyl==0);

CategoriesStruct.diff = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_paired, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 1);


% ######################## COMPARE SAME AND DIFF
goodclasses = [1 2 5]; % those that bring closer
% === the proportion of classes that brings closer vs. further
sylfield = 'same';
Ngoodclass = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
Nbadclass = sum(CategoriesStruct.(sylfield).NAll) - Ngoodclass;
disp(['*** ' sylfield ': ' num2str(Ngoodclass) '/' num2str(Ngoodclass+Nbadclass)]);

sylfield = 'diff';
Ngoodclass = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
Nbadclass = sum(CategoriesStruct.(sylfield).NAll) - Ngoodclass;
disp(['*** ' sylfield ': ' num2str(Ngoodclass) '/' num2str(Ngoodclass+Nbadclass)]);



%% ####################### compare distributions to shuffle distributions
% to shuffle: 1) null hypothesis is that AFP bias is independent of pair
% (i.e. each syl in a pair gets random draw from entire distribution of
% biases (across all datapoints across contexts and syls)
%% SHUFFLE PROCEDURE (COMBINED SAME AND DIFF TO TEST DIFFERENCES BETWEEN THEM)
% ################################################ SAME/DIFF COMBINED
NAll_same = [];
NAll_diff = [];
SepPBS_same = [];
SepMUSC_same = [];

SepPBS_diff = [];
SepMUSC_diff = [];

Ncycles = 2500;
for nn = 1:Ncycles
    
    % --- shuffle all biases (maintaining everything else)
    % NOTE: shuffle at level of syl, not at level of paired syls, since
    % a given syl contributes to multiple pairs
        indshuff = randperm(length(Pitch_MUSC_minus_PBS));
    PitchShifts_shuff = Pitch_MUSC_minus_PBS(indshuff);
    
    % ============================== recalculate Pitch shifts paired
    PitchShiftBothSyls_SHUFF = {};
    NumExpts=max(ExptCounter);
    
    for i=1:NumExpts
        Inds=find(ExptCounter==i);
        % ---- get all pairwise stuff
        for j=1:length(Inds)
            ind1=Inds(j);
            for jj=j+1:length(Inds)
                ind2=Inds(jj);
                % ==== collect pairwise stats
                PitchShiftBothSyls_SHUFF=[PitchShiftBothSyls_SHUFF [PitchShifts_shuff(ind1) PitchShifts_shuff(ind2)]];
            end
        end
    end
    
    % =================================== get distribution of classes
    % ----------- SAME
    inds = IsSameSyl==1;
    structtmp = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_SHUFF, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 0);
    
    NAll_same = [NAll_same; structtmp.NAll'];
    SepPBS_same = [SepPBS_same; structtmp.SeparationPBS'];
    SepMUSC_same = [SepMUSC_same; structtmp.SeparationMUSC'];

    % ----------- DIFF
    inds = IsSameSyl==0;
    structtmp = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_SHUFF, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 0);
    
    NAll_diff = [NAll_diff; structtmp.NAll'];
    SepPBS_diff = [SepPBS_diff; structtmp.SeparationPBS'];
    SepMUSC_diff = [SepMUSC_diff; structtmp.SeparationMUSC'];

   
end
CategoriesStruct.same_SHUFF.NAll_All = NAll_same;
CategoriesStruct.same_SHUFF.SepPBS = SepPBS_same;
CategoriesStruct.same_SHUFF.SepMUSC= SepMUSC_same;

CategoriesStruct.diff_SHUFF.NAll_All = NAll_diff;
CategoriesStruct.diff_SHUFF.SepPBS= SepPBS_diff;
CategoriesStruct.diff_SHUFF.SepMUSC= SepMUSC_diff;


% ############################################ PLOT 
lt_figure; hold on;

% ====== COMPARE EACH POSITION (PROPORTION)
lt_subplot(2,1,1); hold on;

nsame = CategoriesStruct.same.NAll./sum(CategoriesStruct.same.NAll);
ndiff = CategoriesStruct.diff.NAll./sum(CategoriesStruct.diff.NAll);
ndat = nsame - ndiff;

nsame_null = CategoriesStruct.same_SHUFF.NAll_All./sum(CategoriesStruct.same.NAll);
ndiff_null = CategoriesStruct.diff_SHUFF.NAll_All./sum(CategoriesStruct.diff.NAll);
nnull = nsame_null - ndiff_null;

% - plot null
X = 1:length(ndat);
plot(X, nnull, '-r');

% - plot dat
lt_plot_bar(X, ndat);

% - calculate p
for j=1:length(ndat)
    
%     p = sum(nnull(:,j)>ndat(j))/size(nnull,1);
    p = sum(nnull(:,j)>ndat(j))/size(nnull,1);
    lt_plot_text(j-0.3, 1.2*ndat(j), ['p=' num2str(p)], 'b');
    lt_plot_annotation(1, 'one-sided', 'b');
    
end

% ==================== COMPARE "GOOD" CLASSES
nsame = sum(CategoriesStruct.same.NAll(goodclasses)./sum(CategoriesStruct.same.NAll));
ndiff = sum(CategoriesStruct.diff.NAll(goodclasses)./sum(CategoriesStruct.diff.NAll));
ndat = nsame - ndiff;

nsame_null = sum(CategoriesStruct.same_SHUFF.NAll_All(:,goodclasses)./sum(CategoriesStruct.same.NAll),2);
ndiff_null = sum(CategoriesStruct.diff_SHUFF.NAll_All(:,goodclasses)./sum(CategoriesStruct.diff.NAll), 2);
nnull = nsame_null - ndiff_null;

lt_subplot(2,2,3); hold on;
title('diff (good classes');
plot(1, nnull, 'o', 'Color', 'r');
lt_plot_bar(1, ndat);

% -- plot p
    p = sum(nnull>ndat)/size(nnull,1);
    lt_plot_text(1, 1.2*ndat, ['p=' num2str(p)], 'b');
    lt_plot_annotation(1, 'one-sided', 'b');




% ######################################## plot null distribution
lt_figure; hold on;

% ============================================ same
sylfield = 'same';
lt_subplot(3,1,1); hold on;

% ----------- RUN
title(sylfield);
NAll_All = CategoriesStruct.([sylfield '_SHUFF']).NAll_All;

 % --- shuffle
X = 1:length(CategoriesStruct.(sylfield).NAll);
plot(X, NAll_All', '-r');

% -- actual dat
lt_plot_bar(X, CategoriesStruct.(sylfield).NAll, {'Color', [0.6 0.6 0.6]});

% -- get p value for each categorie
for j=1:length(CategoriesStruct.(sylfield).NAll)
    nthis = CategoriesStruct.(sylfield).NAll(j);
    nnull = NAll_All(:,j);
    p = sum(nnull>nthis)/length(nnull);
    lt_plot_text(j-0.3, 1.2*nthis , ['p=' num2str(p)], 'b');
end


% ============================================ same
sylfield = 'diff';
lt_subplot(3,1,2); hold on;

% ----------- RUN
title(sylfield);
NAll_All = CategoriesStruct.([sylfield '_SHUFF']).NAll_All;

 % --- shuffle
X = 1:length(CategoriesStruct.(sylfield).NAll);
plot(X, NAll_All', '-r');

% -- actual dat
lt_plot_bar(X, CategoriesStruct.(sylfield).NAll, {'Color', [0.6 0.6 0.6]});

% -- get p value for each categorie
for j=1:length(CategoriesStruct.(sylfield).NAll)
    nthis = CategoriesStruct.(sylfield).NAll(j);
    nnull = NAll_All(:,j);
    p = sum(nnull>nthis)/length(nnull);
    lt_plot_text(j-0.3, 1.2*nthis , ['p=' num2str(p)], 'b');
end


% ====================================== [proportion in good class]
% [frequency]
lt_subplot(3,2,5); hold on;
xlabel('same -- diff');
ylabel('count');
% --------------------------------------- SAME
sylfield = 'same';
X = 1;

ndat = CategoriesStruct.(sylfield).NAll(goodclasses);
ndat = sum(ndat);
nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');


% --------------------------------------- DIFF
sylfield = 'diff';
X = 2;

ndat = CategoriesStruct.(sylfield).NAll(goodclasses);
ndat = sum(ndat);
nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');

lt_subplot(3,2,5); hold on;
xlabel('same -- diff');
ylabel('count');


% ====================================== [proportion in good class]
% [PROPORTION]
lt_subplot(3,2,6); hold on;
xlabel('same -- diff');
ylabel('fraction');

% ------------------------ SAME
sylfield = 'same';
X = 1;

ndat = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
ndat_tot = sum(CategoriesStruct.(sylfield).NAll);
ndat = ndat/ndat_tot;

nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
nnull_tot =  sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All,2);
nnull = nnull./nnull_tot;

% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');


% ------------------------ SAME
sylfield = 'diff';
X = 2;

ndat = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
ndat_tot = sum(CategoriesStruct.(sylfield).NAll);
ndat = ndat/ndat_tot;

nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
nnull_tot =  sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All,2);
nnull = nnull./nnull_tot;

% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');



%% SHUFFLE PROCEDURE (COMBINED SAME AND DIFF TO TEST DIFFERENCES BETWEEN THEM)
% ################################################ SAME/DIFF COMBINED

% ############################################ PLOT 
lt_figure; hold on;

% ====== COMPARE EACH POSITION (PROPORTION)
lt_subplot(2,1,1); hold on;

nsame = CategoriesStruct.same.NAll./sum(CategoriesStruct.same.NAll);
ndiff = CategoriesStruct.diff.NAll./sum(CategoriesStruct.diff.NAll);
ndat = nsame - ndiff;

nsame_null = CategoriesStruct.same_SHUFF.NAll_All./sum(CategoriesStruct.same.NAll);
ndiff_null = CategoriesStruct.diff_SHUFF.NAll_All./sum(CategoriesStruct.diff.NAll);
nnull = nsame_null - ndiff_null;

% - plot null
X = 1:length(ndat);
plot(X, nnull, '-r');

% - plot dat
lt_plot_bar(X, ndat);

% - calculate p
for j=1:length(ndat)
    
%     p = sum(nnull(:,j)>ndat(j))/size(nnull,1);
    p = sum(nnull(:,j)>ndat(j))/size(nnull,1);
    lt_plot_text(j-0.3, 1.2*ndat(j), ['p=' num2str(p)], 'b');
    lt_plot_annotation(1, 'one-sided', 'b');
    
end

% ==================== COMPARE "GOOD" CLASSES
nsame = sum(CategoriesStruct.same.NAll(goodclasses)./sum(CategoriesStruct.same.NAll));
ndiff = sum(CategoriesStruct.diff.NAll(goodclasses)./sum(CategoriesStruct.diff.NAll));
ndat = nsame - ndiff;

nsame_null = sum(CategoriesStruct.same_SHUFF.NAll_All(:,goodclasses)./sum(CategoriesStruct.same.NAll),2);
ndiff_null = sum(CategoriesStruct.diff_SHUFF.NAll_All(:,goodclasses)./sum(CategoriesStruct.diff.NAll), 2);
nnull = nsame_null - ndiff_null;

lt_subplot(2,2,3); hold on;
title('diff (good classes');
plot(1, nnull, 'o', 'Color', 'r');
lt_plot_bar(1, ndat);

% -- plot p
    p = sum(nnull>ndat)/size(nnull,1);
    lt_plot_text(1, 1.2*ndat, ['p=' num2str(p)], 'b');
    lt_plot_annotation(1, 'one-sided', 'b');



    
    
%% ============= [SHUFFLE] - COMPARE CHANGE IN SEPARATION
lt_figure; hold on;
title('separation');

Y = {};
% ============= actual dat
sepChange_same = cell2mat(CategoriesStruct.same.SeparationMUSC') - ...
    cell2mat(CategoriesStruct.same.SeparationPBS');
Y{1} = sepChange_same;

sepChange_diff = cell2mat(CategoriesStruct.diff.SeparationMUSC') - ...
    cell2mat(CategoriesStruct.diff.SeparationPBS');
Y{2} = sepChange_diff;

ymean = [mean(Y{1}) mean(Y{2})];
ysem = [lt_sem(Y{1}) lt_sem(Y{2})];
% lt_plot_MultDist(Y, [1 2]);
lt_plot([1 2], ymean, {'Errors', ysem, 'Color', 'k'});

% --- overlay null data
numshuffs = size(CategoriesStruct.same_SHUFF.SepPBS,1);
Yshuff = [];
for j=1:numshuffs
sepchange_same_shuff = cell2mat(CategoriesStruct.same_SHUFF.SepMUSC(j,:)) - ...
    cell2mat(CategoriesStruct.same_SHUFF.SepPBS(j,:));

sepchange_diff_shuff = cell2mat(CategoriesStruct.diff_SHUFF.SepMUSC(j,:)) - ...
    cell2mat(CategoriesStruct.diff_SHUFF.SepPBS(j,:));

y = [mean(sepchange_same_shuff) mean(sepchange_diff_shuff)];
plot([1.1 1.9], y, '-', 'Color', [0.7 0.7 0.7]);

Yshuff = [Yshuff; y];
end

% ========== significance
% -- same
p = sum(Yshuff(:,1) > ymean(1))/size(Yshuff,1);
lt_plot_text(1, 1.2*max(Yshuff(:,1)), ['p=' num2str(p)], 'r');

% -- diff
p = sum(Yshuff(:,2) > ymean(2))/size(Yshuff,1);
lt_plot_text(2, 1.2*max(Yshuff(:,2)), ['p=' num2str(p)], 'r');


% -- same minus diff
p = sum((Yshuff(:,1) - Yshuff(:,2)) > (ymean(1) - ymean(2)))/size(Yshuff,1);
lt_plot_text(1.5, 1.2*max(Yshuff(:,2)), ['(vs) p=' num2str(p)], 'r');

    
    
    
%% ################### [SYL1 = starting higher pitch] [change as fraction of separation]
if (0)
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ######################################## same, but using hand label same type
% === plot all pairwise effect (45 degree) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS) [frac of separation]')
ylabel('syl2 shift (MUSC - PBS) [frac of seapration]')

inds=find(IsSameSyl);
XY = cell2mat(PitchShiftBothSyls_paired(inds)');

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));
% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));

% -- convert pitch shift to fraction of separation
Ysep= abs(ffsyl2-ffsyl1);
Ysep = repmat(Ysep', 1, 2);
XY = XY./Ysep;

% --- PLOT
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;
line([-50 50], [-50 50], 'Color','k');

% --- also plot how much greater pitch of syl 2 is (rel syl 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl1 hz - syl2 hz (base)')

Ysep= abs(ffsyl2-ffsyl1);
plot(XY(:,1), Ysep, 'ok');

% ----------------------------
linkaxes(hsplots, 'x');

end
%% OTHER PLOTS RELATIING MUSC EFFECT TO ACOUSTIC DISTANCE

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




%% distributions of changes in separation

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


% ####################################### REPEAT, BUT USING HAND LAB
% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type pairs [HANDLAB] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=IsSameSyl;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')

% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type pairs [HAND LAB] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=IsSameSyl==0;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


%% COMPARE SEPARATION POST VS PRE MUSC

% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

Y=changeInAbsSeparationDuringMUSC_pairs;

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% --- SAME TYPE
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


% --- DIFF TYPE
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


% ################################################### REPLOT, BUT USING
% HAND LABELED
% --- SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type HANDLAB] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=IsSameSyl;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'b');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type HANDLAB] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% --- DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type HANDLAB] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=IsSameSyl==0;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'r');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type HANDLAB] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds=IsSameSyl==0;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'r');


%% ==== directly compare sametype vs. diff type [all datapoints]
Yallall = {}; % cell1 = same; cell2 = diff;

% ============ SAME
inds = IsSameSyl==1;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{1} = X; % separation
Yallall{1} = Y-X; % change in separation
threshold = max(X)+20000;

% ================ DIFF
inds=IsSameSyl==0;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{2} = X(X<=threshold);
Yallall{2} = Y(X<=threshold) - X(X<=threshold);

% ========================================== PLOIT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

lt_plot_MultDist(Yallall, [1 2]);

%% ==== directly compare sametype vs. diff type, [constraining to those that start similar]
Yallall = {}; % cell1 = same; cell2 = diff;

% ============ SAME
inds = IsSameSyl==1;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{1} = X; % separation
Yallall{1} = Y-X; % change in separation
threshold = max(X);

% ================ DIFF
inds=IsSameSyl==0;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{2} = X(X<=threshold);
Yallall{2} = Y(X<=threshold) - X(X<=threshold);

% ========================================== PLOIT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

lt_plot_MultDist(Yallall, [1 2]);


%% %%%%%%%%%%%%%%%%%%%%%%%% change in separation as a function of Motif position
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ====================== SAME MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same motif');
Yall ={};

% ----- same pair
inds = IsSameMotif==1 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% ----- diff pair
inds = IsSameMotif==1 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2])



% ====================== DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff motif');
Yall ={};

% ----- same motif, same pair
inds = IsSameMotif==0 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% ----- same motif, diff pair
inds = IsSameMotif==0 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2])


%% %%%%%%%%%%%%%%%%%%%%%%%% change in separation as a function of Motif position
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ====================== SAME SYL
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same SYL');
xlabel('same motif -- diff motif');
ylabel('change (ff');
Yall ={};

% ----- same motif
inds = IsSameMotif==1 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% ----- diff motif
inds = IsSameMotif==0 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2])

ylim([-120 120]);
lt_plot_zeroline;

% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(1.5, 0.8, ['(vs)p=' num2str(p)], 'b')



% ====================== DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff syl');
Yall ={};

% -----
inds = IsSameMotif==1 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% -----
inds = IsSameMotif==0 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2]);

ylim([-120 120]);
lt_plot_zeroline;

% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(1.5, 0.8, ['(vs)p=' num2str(p)], 'b')




%% ####################################### SCATTER PLOT [SEPARATION]
ChangeInSep = abs(separationMUSC_paired) - abs(separationPBS_paired);
% ======== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==1;

x = abs(separationPBS_paired(inds));
y = ChangeInSep(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
ylim([-100 100]);
aoctool(x, y, samemotif)

% ======== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff (bk = samemotif; rd = diff)');
inds = IsSameSyl==0;

x = abs(separationPBS_paired(inds));
y = ChangeInSep(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');
lt_plot_zeroline;
ylim([-100 100]);



%% ===== same diff motif (45 degree scatter) [SEPARATION]

% ####################################### SCATTER PLOT
ChangeInSep = abs(separationMUSC_paired) - abs(separationPBS_paired);
% ======== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==1;

x = abs(separationPBS_paired(inds));
y = abs(separationMUSC_paired(inds));
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;

% ylim([-100 100]);

aoctool(x, y, samemotif)



% ======== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==0;

x = abs(separationPBS_paired(inds));
y = abs(separationMUSC_paired(inds));
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
% ylim([-100 100]);

aoctool(x, y, samemotif)


%% ######################33 [SEAPRATION], DEPENDENCE ON  POSOTIION WITHIN MOTIF
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


ChangeInSep = abs(separationMUSC_paired) - abs(separationPBS_paired);

inds = IsSameMotif==1 & ~isnan(NumSylsInBetween);

nsyls = NumSylsInBetween(inds);
syltype = IsSameSyl(inds);
Ycorrchange = ChangeInSep(inds);

% ========================== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
ylabel('separationFF (MUSC minus PBS');
xlabel('num syls in between');
x = nsyls(syltype==1);
y = Ycorrchange(syltype==1);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);
mean(y(x==0)), mean(y(x~=0));

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);

% ============================= DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');

x = nsyls(syltype==0);
y = Ycorrchange(syltype==0);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);

%% ===== MUSC drive same-type to be more correlated? [song corr]
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% =========== SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

% ----- PLOT 1- SCATTER
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);

% ---- PLOT 2 - PAIRED TEST
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');

x = [1 2];
y = [CorrSong_PBS_pairs(inds)' CorrSong_MUSC_pairs(inds)'];
plot(x, y, '-ok');
plot(2, y(:,2), 'or');
% --- plot means
ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.1, ymean, {'Errors', ysem, 'LineStyle', '-', 'Color', 'y'});
% --- sign rank test
p = signrank(y(:,1), y(:,2));
lt_plot_pvalue(p, 'srank', 1);
xlim([0 3]);
ylim([-0.5 1]);
lt_plot_zeroline;


% ======================== DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

% ---------- PLOT 1 - SCATTER
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);

% -------- PLOT 2 - PAIRED TEST
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');

x = [1 2];
y = [CorrSong_PBS_pairs(inds)' CorrSong_MUSC_pairs(inds)'];
plot(x, y, '-ok');
plot(2, y(:,2), 'or');
% --- plot means
ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.1, ymean, {'Errors', ysem, 'LineStyle', '-', 'Color', 'y'});
% --- sign rank test
p = signrank(y(:,1), y(:,2));
lt_plot_pvalue(p, 'srank', 1);
xlim([0 3]);
ylim([-0.5 1]);
lt_plot_zeroline;



% =========================== COMBINED SCATTER PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME(bu), DIFF(rd)');
xlabel('PBS - paired corr [song]');
ylabel('MUSC - paired corr [song');

% 1) SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);


% 2) DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds);
lt_plot_45degScatter(X, Y, color);
xlim([-1 1]); ylim([-1 1]);


% =========================== CHANGE IN CORR RELATIVE TO STARTING CORR
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME(bu), DIFF(rd)');
xlabel('PBS - paired corr [song]');
ylabel('CHANGE IN CORR (musc - pbs)');

Xallall = [];
Yallall = [];
Grpall = [];

% 1) SAME TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_regress(Y, X, 1, 0, 1, 1, color, 0);
lt_plot_zeroline;

Xallall = [Xallall X];
Yallall = [Yallall Y];
Grpall = [Grpall 1*ones(1,length(X))];


% 1) DIFF TYPE
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_PBS_pairs(inds);
Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_regress(Y, X, 1, 0, 1, 1, color, 0);

Xallall = [Xallall X];
Yallall = [Yallall Y];
Grpall = [Grpall 2*ones(1,length(X))];


% ######################################### ANCOVA (change in intercept)
[h,atab,ctab,stats] = aoctool(Xallall', Yallall', Grpall', 0.05, ...
    '','','','on','parallel lines');
c = multcompare(stats)



%% ############################# difference distyributions
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('corr(MUSC)-corr(PBS)');
X=-1:0.05:1;

Yallall ={};
% same
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

Yallall{1} = Y;

% diff
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

Y=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
lt_plot_histogram(Y, X, 1, 1, '', 1, color);

Yallall{2} = Y;



% ================== Plot distributions side by side
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('change in corr (song)');

lt_plot_MultDist(Yallall, [1 2])
lt_plot_zeroline;

%% ###################### CORR, dependence on WHETHER same motif)


% ================================ SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==1;

x = CorrSong_PBS_pairs(inds);
y = CorrSong_MUSC_pairs(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
line([-0.5 1], [-0.5 1])
% ylim([-100 100]);

aoctool(x, y, samemotif)


% ======== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff (bk = samemotif; rd = diff)');
xlabel('PBS corr');
ylabel('MUSC corr');
inds = IsSameSyl==0;

x = CorrSong_PBS_pairs(inds);
y = CorrSong_MUSC_pairs(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
line([-0.5 1], [-0.5 1])
% ylim([-100 100]);

aoctool(x, y, samemotif)


%% ======= COMPARE SAME AND DIFF ... (CORR)
x = CorrSong_PBS_pairs;
y = CorrSong_MUSC_pairs;
samesyl = IsSameSyl;

aoctool(x, y, samesyl)


%% ===================== CORR, SEPARATE BY MOTIF (DIRECTL YCOMPARE SASME AND DIFF)

ChangeInCorr = CorrSong_MUSC_pairs - CorrSong_PBS_pairs;

% ===================================== SAME MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME MOTIF (bk = same syl; rd = diff syl)');
Yall = {};

% ---- same syl
inds = IsSameMotif==1 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff syl
inds = IsSameMotif==1 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)


% ===================================== DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF MOTIF (bk = same syl; rd = diff syl)');
Yall = {};

% ---- same syl
inds = IsSameMotif==0 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff syl
inds = IsSameMotif==0 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)


%% ====================== CORR, COMPARE WITHIN OR ACROSS MOTIFS
ChangeInCorr = CorrSong_MUSC_pairs - CorrSong_PBS_pairs;

% ========== SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL (left = same motif; rt = diff motif)');
Yall = {};
ylabel('change in corr (MUSC-PBS)');

% ---- same motif
inds = IsSameMotif==1 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff motif
inds = IsSameMotif==0 & IsSameSyl==1;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)
lt_plot_zeroline
ylim([-1 1]);
% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i-1, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(0.5, 0.8, ['(vs)p=' num2str(p)], 'b')

% ========== DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL (bk = same motif; rd = diff motif)');
Yall = {};

% ---- same motif
inds = IsSameMotif==1 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{1} = y;

% ---- diff motif
inds = IsSameMotif==0 & IsSameSyl==0;
y = ChangeInCorr(inds);
Yall{2} = y;

lt_plot_MultDist(Yall, [0 1], 1)
lt_plot_zeroline
ylim([-1 1]);
% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i-1, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(0.5, 0.8, ['(vs)p=' num2str(p)], 'b')

%% ######################33 CORR, DEPENDENCE ON MOTIF POSOTIION
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


ChangeInCorr = CorrSong_MUSC_pairs - CorrSong_PBS_pairs;

inds = IsSameMotif==1 & ~isnan(NumSylsInBetween);

nsyls = NumSylsInBetween(inds);
syltype = IsSameSyl(inds);
Ycorrchange = ChangeInCorr(inds);

% ========================== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
ylabel('corr (MUSC minus PBS');
xlabel('num syls in between');
x = nsyls(syltype==1);
y = Ycorrchange(syltype==1);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);
mean(y(x==0)), mean(y(x~=0));

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);

% ============================= DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');

x = nsyls(syltype==0);
y = Ycorrchange(syltype==0);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
ylim([-1 1]);

%% ===== MUSC drive same-type to be more correlated? [motif corr]
if (0) % REASON: focus on song corr above - can have more datapoints
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
end


%% ==== effect of MUSC on corr correlate with effect of MUSC on separation?
% --- SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('corr(MUSC) - corr(PBS)[song]');
ylabel('change in abs separation');
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
color='b';

X=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

% --- diff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('corr(MUSC) - corr(PBS)[song]');
ylabel('change in abs separation');
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
color='r';

X=CorrSong_MUSC_pairs(inds)-CorrSong_PBS_pairs(inds);
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)


disp('Effect of MUSC on pairwise corr did not predict effect of MUSC on separation [both sim and diff]');


%% === change in corr versus corr
if (0) % THIS IS REPLACED BY ABOVE (correaltion of change in corr vs. starting corr)
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
end

%% ==== distributions of correlations (MUSC and PBS);

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

% -- SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end
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


% -- DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
xlabel('corr(song)');
lt_plot_annotation(1, 'dash=PBS; sol:MUSC', 'k');
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end
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


% === ALL PAIRS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
inds=find(AcousticDist_pairs>0);

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


% === SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same]');
color='b';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff]');
color='r';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


%% ==== change in acoustic separation (vs. MUSC acoustic sep)
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === ALL
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

% === SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same]');
color='b';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end

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
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

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
assert(length(inda)>length(tvals1)/2, 'why not much intersect?')

corrsong = corr(ffvals1', ffvals2');
assert(~any(isnan(corrsong(:))));
end

function structout = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_paired, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, ...
    plotON)

NAll = [];
NlessSepAll = [];
PitchShiftsAll = {};
PitchPBSAll = {};
SeparationPBS = {};
SeparationMUSC = {};

% -------------------- EXTRACT THINGS
XY = cell2mat(PitchShiftBothSyls_paired(inds)');
YsepPBS = abs(separationPBS_paired(inds));
if (0)
YsepMUSC = abs(separationMUSC_paired(inds));
else
YsepMUSC = Pitch_PBSraw(Pairs_OriginalInds(inds, :)) + XY;
YsepMUSC = abs(diff(YsepMUSC,1,2))';
end

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));

% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));
ffsyl2COPY = ffsyl2;
ffsyl2(indflip) = ffsyl1(indflip);
ffsyl1(indflip) = ffsyl2COPY(indflip);
clear ffsyl2COPY


% -- annotate cases where MUSC brought separation closer
indsMuscCloser = YsepMUSC < YsepPBS;

% ================ GO THRU CATEGORIES ONE BY ONE AND EXTRACT 1) N, 2) HOW
% MANY DECREASED SEPARATION, 3) PITCH SHIFTS
% ########## 1) down, up [i.e. syl 1 down, syl 2 up]
indtmp = XY(:,1) < 0 & XY(:,2) > 0 ;
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 2) up(less), up(more)
indtmp = XY(:,1) > 0 & XY(:,2) > 0 & XY(:,1)<XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 2) up(more), up(less)
indtmp = XY(:,1) > 0 & XY(:,2) > 0 & XY(:,1)>XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 3) dn(less), dn(more)
indtmp = XY(:,1) < 0 & XY(:,2) < 0 & XY(:,1)>XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 3) dn(more), dn(less)
indtmp = XY(:,1) < 0 & XY(:,2) < 0 & XY(:,1)<XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 4) up, dn
indtmp = XY(:,1) > 0 & XY(:,2) < 0 ;
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ################################### PLOT
if plotON==1
lt_figure; hold on;

% =========== 1) bar plot of N(with fill showing how many reduced
% separation)
lt_subplot(4,1,1); hold on;
ylabel('Num cases (red = sep decrease)');
X = 1:length(NAll);
lt_plot_bar(X, NAll, {'Color', [0.6 0.6 0.6]});
lt_plot_bar(X, NlessSepAll, {'Color', 'r'});

% =========== 2) pitch shifts
lt_subplot(4,1,2); hold on
ylabel('pitch shift (hz');
for j=1:length(PitchShiftsAll)
   
    plot([j-0.2 j+0.2], PitchShiftsAll{j}','-ok');

end
lt_plot_zeroline

% =========== 3) pitch (PBS) -- pitch (MUSC)
lt_subplot(4,1,3); hold on
% ylabel('pitch shift (hz');
for j=1:length(PitchShiftsAll)
    
    numcases = size(PitchShiftsAll{j},1);
    plotcols = lt_make_plot_colors(numcases, 0, 0);
    for jj=1:numcases
    
    ffPBS = PitchPBSAll{j}(jj,:)';
    ffPBS = ffPBS-mean(ffPBS); % subtract mean
    ffMUSC = ffPBS + PitchShiftsAll{j}(jj,:)';
    
    % plot 2 contexts separate x places
%      plot([j-0.2 j+0.2], [ffPBS ffMUSC]', '-', 'Color', plotcols{jj});
    plot([j-0.3 j-0.2], [ffPBS(1) ffMUSC(1)]', '-', 'Color', plotcols{jj}, 'LineWidth', 2);
    plot([j+0.2 j+0.3], [ffPBS(2) ffMUSC(2)]', '-', 'Color', plotcols{jj}, 'LineWidth', 2);
    plot([j-0.2 j+0.2], [ffMUSC(1) ffPBS(2)], ':', 'Color', [0.7 0.7 0.7]);
        
    end
    
%     
%     plot([j-0.2 j+0.2], PitchShiftsAll{j}','-ok');

end
lt_plot_zeroline


% =========== 4) separtiaon (PBS) -- sep (MUSC)
lt_subplot(4,1,4); hold on
title('seprataiopn');
for j=1:length(PitchShiftsAll)
   
    yPBS = SeparationPBS{j};
    yMUSC = SeparationMUSC{j};
    plot([j-0.2 j+0.2], [yPBS' yMUSC']','-ok');

end
lt_plot_zeroline
end


% ================== OUTPUT
structout.NAll = NAll;
structout.NlessSepAll = NlessSepAll;
structout.PitchShiftsAll = PitchShiftsAll;
structout.PitchPBSAll = PitchPBSAll;
structout.SeparationPBS =SeparationPBS;
structout.SeparationMUSC = SeparationMUSC;


end