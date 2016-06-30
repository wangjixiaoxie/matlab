function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANpreprocess(SeqDepPitch_AcrossBirds, PARAMS, musc_day_window_WNdayinds, musc_day_window_Bidirdayinds_Bidir, pu53_use_later_days, debugON, single_targ_UseMaintainedEarly)
% musc_day_window_WNdayinds=[3 10]; % from WN day 3 to 10, collect all for analysis
% pu53_use_later_days=1; then for (birdname, 'pu53wh88') & strcmp(exptname, 'SeqDepPitchLMAN'), uses later days, because inactivation did not work for earlier days


%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

%% 2) WHAT ARE DAYS WITH MUSCMOL?

for i=1:NumBirds;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        muscdays_inds=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow));
        
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds=muscdays_inds;
        
    end
end





%% 3) EXTRACT ACROSS DAY MEANS OF EFFECT OF MUSC [SINGLE DIR]
disp(' -- SINGLE DIR EXTRACTION ');
LagTimes_All=[];
PBSstartTimeALL=[];
StartTimeRangeALL=[];
StartTimeMeanALL=[];
StartTimeSTDALL=[];

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        muscdays_inds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        
        WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        
        % --- start of window in inds from start of expt
        WindowFirstDay=WNonInd+musc_day_window_WNdayinds(1)-1;
        WindowLastDay=WNonInd+musc_day_window_WNdayinds(2)-1;
        
        % --- go through all days in the window that are also musc days
        % - extract raw data
        MuscDaysWithinWindow=muscdays_inds(muscdays_inds>=WindowFirstDay & muscdays_inds<=WindowLastDay);
        
        % ------ Only keep single dir days
        daybeforemulti=1000; % 1000 will always be larger than window
        daybeforeSame=1000;
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            daybeforemulti=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        end
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart_Ind')
            daybeforeSame=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind;
        end
        
        LastGoodDay=min([daybeforemulti daybeforeSame]);
        
        if strcmp(birdname, 'rd28pu64') & strcmp(exptname, 'SeqDepPitchLMAN2');
            % then actually did samedir early, but not annotated (becuase
            % only started hitting all variants late)
            LastGoodDay=20;
        end
        
        
        MuscDaysWithinWindow=MuscDaysWithinWindow(MuscDaysWithinWindow<=LastGoodDay);
        
        
        % ----------
        if strcmp(birdname, 'pu53wh88') & strcmp(exptname, 'SeqDepPitchLMAN');
            if pu53_use_later_days==1;
                % then for this expt, use later days
                MuscDaysWithinWindow=intersect(muscdays_inds, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.days_consolid_early);
            elseif pu53_use_later_days==2
                % then throw out this experiment (as no inactivation until
                % late days)
                continue
            end
        end
        
        if debugON==1
            lt_figure; hold on;
            title([birdname '-' exptname ' SING DIR (with sametype)']);
            xlim([0 40]);
            line(xlim, [0 0]);
        end
        
        
        % ==== IF WANT TO CONSTRAIN relative to maintained shift:
        % window = [min day (as above), min(max_day_as_above, day 4 of maintained
        % shift)]
        if single_targ_UseMaintainedEarly==1
            % find consolidation day for this expt
            ind1=find(strcmp(PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, birdname));
            ind2=find(strcmp(PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, exptname));
            ind3=intersect(ind1+2, ind2+1);
            if ~isempty(ind3)
                tmp=PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind3};
                maintained_shift_day3=tmp(1)+3;
                maintained_shift_day3=maintained_shift_day3+WNonInd-1;
                disp([birdname '-' exptname ': OLD MUSC DAYS: ' num2str(MuscDaysWithinWindow) '; NEW MAX DAY: ' num2str(maintained_shift_day3)]);
                MuscDaysWithinWindow=MuscDaysWithinWindow(MuscDaysWithinWindow<=maintained_shift_day3);
            end
        end
        
        % ---- IF NO DAYS, THEN HAVE MESASGE
        if isempty(MuscDaysWithinWindow)
            disp(['WARNING WARNING: no musc days found in window for ' birdname '-' exptname '- vals will be empty/nan']);
            
        end
        
        % === lag time collect
        lagtime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.Lag_time;
LagTimes_All=[LagTimes_All lagtime];

% === start of PBS collection time
pretime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.PBS_window(1);
medianSwitchTime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_MedianStartTime;
PBSstartTime=medianSwitchTime+pretime;
PBSstartTimeALL=[PBSstartTimeALL PBSstartTime];

% ==== start of MUSC switch time
StartTimes=[];
for k=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds);
    if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds{k})
        continue
    end
    stime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds{k};
StartTimes=[StartTimes stime.start];
end
StartTimeRange=max(StartTimes)-min(StartTimes);
StartTimeMean=mean(StartTimes);
StartTimeSTD=std(StartTimes);

StartTimeRangeALL=[StartTimeRangeALL StartTimeRange];
StartTimeMeanALL=[StartTimeMeanALL StartTimeMean];
StartTimeSTDALL=[StartTimeSTDALL StartTimeSTD];


        
        % --- EXTRACT DATA FOR THOSE DAYS
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            
            FFpbsAll=[];
            TvalsPBSall=[];
            FFmuscAll=[];
            TvalsMUSCall=[];
            if ~isempty(MuscDaysWithinWindow)
                for dayInd=MuscDaysWithinWindow
                    
                    ffpbs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{dayInd};
                    ffmusc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{dayInd};
                    
                    FFpbsAll=[FFpbsAll ffpbs];
                    FFmuscAll=[FFmuscAll ffmusc];
                    
                    TvalsPBSall=[TvalsPBSall SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{dayInd}];
                    TvalsMUSCall=[TvalsMUSCall SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{dayInd}];
                    
                    if debugON==1
                        
                        if strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl, syl);
                            
                            plot(dayInd, mean(FFpbsAll), 'ok'); % then is first targ
                            plot(dayInd, mean(FFmuscAll), 'sk');
                        elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ
                            
                            plot(dayInd, mean(FFpbsAll), 'ob');
                            plot(dayInd, mean(FFmuscAll), 'sb');
                        end
                    end
                end
                
            end
            
            
            
            
            
            
            % ---- Put into output for this syl
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds=MuscDaysWithinWindow;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).FFminusBaseWithinWindow_PBS_vals=FFpbsAll;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_pbs=mean(FFpbsAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_pbs=lt_sem(FFpbsAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).TvalsWithinWindow_PBS=TvalsPBSall;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).FFminusBaseWithinWindow_MUSC_vals=FFmuscAll;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_musc=mean(FFmuscAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_musc=lt_sem(FFmuscAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).TvalsWithinWindow_MUSC=TvalsMUSCall;
        end
        
        % ===== OUTPUT to params
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.final_extracted_window.Window_WNdayinds=musc_day_window_WNdayinds;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.final_extracted_window.dayIndsUsed=MuscDaysWithinWindow;
        
        
        disp([birdname '-' exptname ' used dayInds: ' num2str(MuscDaysWithinWindow) '; WNonInd: ' num2str(WNonInd)]);
        
    end
end



% ======= PLOT lag time and start time
lt_figure; 
lt_plot_annotation(1, ['lag time(std) = ' num2str(mean(LagTimes_All)) '(' num2str(std(LagTimes_All)) ')'], 'k');
lt_plot_annotation(2, ['PBS start time(std) = ' num2str(mean(PBSstartTimeALL)) '(' num2str(std(PBSstartTimeALL)) ')'], 'b');
lt_plot_annotation(3, ['MUSC start time(std) = ' num2str(mean(StartTimeMeanALL)) '(' num2str(std(StartTimeMeanALL)) ')'], 'b');
lt_plot_annotation(4, ['MUSC start time STD (mean and largest) = ' num2str(mean(StartTimeSTDALL)) '(' num2str(max(StartTimeSTDALL)) ')'], 'r');


%% === extract for bidir

disp(' --');
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % ==== continue if has bidir
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            %             disp([birdname '-' exptname]);
            continue;
        end
        
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        muscdays_inds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        
        MultidirDayOneInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
        
        % --- start of window in inds from start of multidir
        WindowFirstDay=MultidirDayOneInd+musc_day_window_Bidirdayinds_Bidir(1)-1;
        WindowLastDay=MultidirDayOneInd+musc_day_window_Bidirdayinds_Bidir(2)-1;
        
        
        % --- go through all days in the window that are also musc days
        % - extract raw data
        MuscDaysWithinWindow=muscdays_inds(muscdays_inds>=WindowFirstDay & muscdays_inds<=WindowLastDay);
        %         disp(WindowFirstDay);
        %         disp(['original musc days: ' num2str(muscdays_inds)]);
        % ---- throw out days after WN end
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            
            WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            %             disp(WNoffInd)
            MuscDaysWithinWindow=MuscDaysWithinWindow(MuscDaysWithinWindow<=WNoffInd);
        end
        
        % ---- IF NO DAYS, THEN HAVE MESASGE
        if isempty(MuscDaysWithinWindow)
            disp(['WARNING: no BIDIR musc days found in window for ' birdname '-' exptname '- vals will be empty/nan']);
            %         else
            %             disp(['Bidir days (overal inds) -' birdname '-' exptname ': ' num2str(MuscDaysWithinWindow)]);
        end
        
        
        
        
        
        % --- EXTRACT DATA FOR THOSE DAYS
        if debugON==1
            
            lt_figure; hold on;
            title([birdname '-' exptname 'BIDIR: (note: is cumulative pitch over days']);
            xlim([0 50]);
            line(xlim, [0 0]);
        end
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            
            FFpbsAll=[];
            TvalsPBSall=[];
            FFmuscAll=[];
            TvalsMUSCall=[];
            if ~isempty(MuscDaysWithinWindow)
                for dayInd=MuscDaysWithinWindow
                    
                    FFpbsAll=[FFpbsAll SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{dayInd}];
                    FFmuscAll=[FFmuscAll SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{dayInd}];
                    
                    TvalsPBSall=[TvalsPBSall SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{dayInd}];
                    TvalsMUSCall=[TvalsMUSCall SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{dayInd}];
                    
                    
                    if debugON==1
                        
                        if any(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls, syl));
                            if strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl, syl);
                                
                                plot(dayInd, mean(FFpbsAll), 'ok'); % then is first targ
                                plot(dayInd, mean(FFmuscAll), 'sk');
                            else
                                plot(dayInd, mean(FFpbsAll), 'ob');
                                plot(dayInd, mean(FFmuscAll), 'sb');
                            end
                        end
                    end
                end
            end
            
            
            
            
            % ---- Put into output for this syl
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.dayInds=MuscDaysWithinWindow;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).FFminusBaseWithinWindow_PBS_vals=FFpbsAll;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).meanFF_pbs=mean(FFpbsAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).semFF_pbs=lt_sem(FFpbsAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).TvalsWithinWindow_PBS=TvalsPBSall;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).FFminusBaseWithinWindow_MUSC_vals=FFmuscAll;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).meanFF_musc=mean(FFmuscAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).semFF_musc=lt_sem(FFmuscAll);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window_bidir.(syl).TvalsWithinWindow_MUSC=TvalsMUSCall;
        end
        
        
        % ===== OUTPUT to params
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.final_extracted_window_bidir.Window_WNdayinds=musc_day_window_WNdayinds;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.final_extracted_window_bidir.dayIndsUsed=MuscDaysWithinWindow;
        
        
        disp(['BIDIR: ' birdname '-' exptname ' used dayInds: ' num2str(MuscDaysWithinWindow)]);
        
        
        
    end
end


%         % --- EXTRACT DATA FOR THOSE DAYS
%         SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         for j=1:length(SylsUnique)
%             syl=SylsUnique{j};
%
%             FFpbsAll=[];
%             TvalsPBSall=[];
%             FFmuscAll=[];
%             TvalsMUSCall=[];
%             if ~isempty(MuscDaysWithinWindow)
%                 for dayInd=MuscDaysWithinWindow
%
%                     FFpbsAll=[FFpbsAll SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{dayInd}];
%                     FFmuscAll=[FFmuscAll SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{dayInd}];
%
%                     TvalsPBSall=[TvalsPBSall SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{dayInd}];
%                     TvalsMUSCall=[TvalsMUSCall SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{dayInd}];
%
%                 end
%             end
%
%
%
%
%             % ---- Put into output for this syl
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds=MuscDaysWithinWindow;
%
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).FFminusBaseWithinWindow_PBS_vals=FFpbsAll;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_pbs=mean(FFpbsAll);
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_pbs=lt_sem(FFpbsAll);
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).TvalsWithinWindow_PBS=TvalsPBSall;
%
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).FFminusBaseWithinWindow_MUSC_vals=FFmuscAll;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_musc=mean(FFmuscAll);
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_musc=lt_sem(FFmuscAll);
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).TvalsWithinWindow_MUSC=TvalsMUSCall;
%         end
%
%         % ===== OUTPUT to params
%         SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.final_extracted_window.Window_WNdayinds=musc_day_window_WNdayinds;
%         SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.final_extracted_window.dayIndsUsed=MuscDaysWithinWindow;
%
%
%         disp([birdname '-' exptname ' used dayInds: ' num2str(MuscDaysWithinWindow) '; WNonInd: ' num2str(WNonInd)]);
%
%     end
% end


