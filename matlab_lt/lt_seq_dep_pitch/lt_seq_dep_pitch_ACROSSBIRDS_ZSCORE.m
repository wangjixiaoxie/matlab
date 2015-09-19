function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_ZSCORE(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 8/13/15 - Analyzes z-scored learning data.
% e.g. params:
% PARAMS.zscore.zthresh_learning=1; % take 1st N days after target passes this value of shift (will take max of default day (see below) or this day)
% PARAMS.zscore.default_min_day=3; 
% PARAMS.zscore.max_day=3; % inclusive, limit days to within min and max days (wn start =1) then throw out this experiment
% PARAMS.zscore.N_days_post_thresh=1; % how many days to use to quantify learning
% Params.zscore.account_for_NoData=1; % If early WN days don't have data, then the window will be moved forward
% % the number of no-data days (assumes no singing early days)

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% EXTRACT Z-SCORED LEARNING METRICS

zthresh_learning=PARAMS.zscore.zthresh_learning; % take 1st N days after target passes this value of shift (will take max of default day (see below) or this day)
default_min_day=PARAMS.zscore.default_min_day; 
max_day=PARAMS.zscore.max_day; % if it takes this num days or longer, then throw out this experiment
N_days_post_thresh=PARAMS.zscore.N_days_post_thresh; % how many days to use to quantify learning

% ===== FOR EACH EXPERIMENT FIND DAYS TO USE
lt_figure; hold on;
title('Each expt, targ syl, zscore learning (mean across days pass thr) and day num (rel to WN start)');
ylabel('z-score learning');
xlabel('day num');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % OVERWRITE OLD DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}, 'Data_ZSCORE');
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}=rmfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}, 'Data_ZSCORE');
        end
        
        % ===== TARGET: WHEN DOES IT PASS CRITERIA FOR LEARNING?
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        Zvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_zscore; % learning trajectory of target
        
        % 1st day pass threshold (limiting search to days within the min
        % and max days)
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        MinDay=WNday1+default_min_day-1;
        MaxDay=WNday1+max_day-1; % convert from WN day to all day inds
        
        % === 1st WN day no data? if so how many consecutive days
        % (including wn day 1) without data? - add that to WN window
        if PARAMS.zscore.account_for_NoData==1;
            NoDataDays=find(cellfun(@isempty, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals));
            
            % how many continuous days are empty starting and including 1st
            % WN day?
            if any(NoDataDays==WNday1);
                num_empty_days=1;
                
                ind=find(NoDataDays==WNday1); % corresponds to 1st WN day
                
                for j=ind+1:length(NoDataDays);
                    
                    if NoDataDays(j)==WNday1+j-ind;
                        num_empty_days=j-ind+1;
                    else
                        % a day has data - stop
                        break
                    end
                    
                end
                
                disp(['ADDED ' num2str(num_empty_days) ' (i.e. missing starting WN days) to ' birdname '-' exptname]);
                MinDay=MinDay+num_empty_days;
                MaxDay=MaxDay+num_empty_days;
            end
        end
                
        % ===
        % make sure WN duration for this experiment was long enough to even
        % have data within the desired day window
        if length(Zvals)<MaxDay;
            disp('PROBLEM - not enough WN days to even contain the day window you asked for');
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds=nan;
            continue;
        end
        
            
            
            
        FirstDay_PassThresh=find(abs(Zvals(MinDay:MaxDay))>zthresh_learning, 1)+MinDay-1; % over all days
        
        % If is empty, then no day in this window passes threshold. don't
        % give a value;
        if isempty(FirstDay_PassThresh);
                    % then is problem
            disp(['PROBLEM - ' birdname '-' exptname ' does not pass learning threshold within days window']);
            FirstDay_PassThresh=nan;
    
        end
        
        
%         FirstDay_PassThresh=find(abs(Zvals)>zthresh_learning, 1);
%         
        % If is before default min day, then default to min day
%         WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
%         WNday3=WNday1+default_min_day-1;
%         
%         if FirstDay_PassThresh<WNday3;
%             FirstDay_PassThresh=WNday3;
%         end
        
        % make sure does not go past WN days
        WNday_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        
        if FirstDay_PassThresh+N_days_post_thresh-1>WNday_last;
            % then is problem, not enough days
            disp(['PROBLEM - ' birdname '-' exptname ' passes thr after WN is over']);
            FirstDay_PassThresh=nan;
        end
        
%         % make sure is not past max day (i.e. takes too long)
%         max_day_ind=WNday1+max_day-1; % convert from WN day to all day inds
%         if FirstDay_PassThresh>=max_day_ind;
%             disp(['PROBLEM - ' birdname '-' exptname ' passes z-score thresh too late']);
%             
%             FirstDay_PassThresh=nan;
%         end            

        % make sure enough days overall
        if FirstDay_PassThresh+N_days_post_thresh-1>length(Zvals);
            disp(['PROBLEM - ' birdname '-' exptname ' simply does not have enough days post pass-threshold']);
            
            FirstDay_PassThresh=nan;
        end   
        
        % make sure no day has nans
        if ~isnan(FirstDay_PassThresh)
        for j=FirstDay_PassThresh:FirstDay_PassThresh+N_days_post_thresh-1;
            if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_zscore(j));
                disp(['PROBLEM - ' birdname '-' exptname ' one day in learning window has nan']);
               FirstDay_PassThresh=nan;
            end
        end
        end
        
        
        % ======== OUTPUT DATA TO STRUCTURE
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds=FirstDay_PassThresh:FirstDay_PassThresh+N_days_post_thresh-1;
        
        % ------ GO TO NEXT EXPERIMENT IF FIRST DAY IS NAN - i.e. doesn't
        % pass criterion
        if isnan(FirstDay_PassThresh);
            continue;
        end
        
                
        % --- PLOT: scatter of days used (an corresponding z-scores) to
        % quantify learning
        FirstDay_PassThresh_WNdays=FirstDay_PassThresh-WNday1+1;
        lt_plot(FirstDay_PassThresh_WNdays, mean(Zvals(FirstDay_PassThresh:FirstDay_PassThresh+N_days_post_thresh-1)), {'Color', 'r'});
        lt_plot_text(FirstDay_PassThresh_WNdays+0.1, mean(Zvals(FirstDay_PassThresh:FirstDay_PassThresh+N_days_post_thresh-1)), [birdname '-' exptname]);
        
        Xlim=xlim;
        xlim([0 max_day+2]);
        
        lt_plot_zeroline;
        plot(xlim, [zthresh_learning zthresh_learning]);
        plot(xlim, [-zthresh_learning -zthresh_learning]);
        
        
        % ====== FOR ALL SYLS, EXTRACT LEARNING METRIC USING THOSE DAYS
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Days=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            ffvals=[];
            tvals=[];
            for k=Days;
            % ----- Collect all ffvals and convert each to a z-score
            ffvals=[ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
            tvals=[tvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{k})];
            end
            
            % --- convert to z-score rel to baseline.
            baseline_meanFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            baseline_stdFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
            
            Zvals=(ffvals-baseline_meanFF)./baseline_stdFF;
            
            % ========= OUTPUT DATA - Save z-score values
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).FFvals=ffvals;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).Zvals=Zvals;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).mean_Zscore=mean(Zvals);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).sem_Zscore=lt_sem(Zvals);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).Tvals=tvals;
            
        end
    end
end

disp('NOTE: if any PROBLEM came, up the experiment was given a nan for FirstDay_PassThresh');
       


%% Plot learning across days - zscored

count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);     
        title([birdname '-' exptname]);
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        plotcols=lt_make_plot_colors(length(syls_unique), 0, 0);
        hplot=[];
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            
            
            % ================ GATHER Z-SCORE LEARNING
            Zvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore;
            
            
            % +++++++++++++++++++
            % ==== PLOT
            hplot(j)=plot(1:length(Zvals), Zvals, 'Color',plotcols{j});
            lt_plot(1:length(Zvals), Zvals, {'Color',plotcols{j}});
        end
        
        legend(hplot, syls_unique);
        line(xlim, [1.5 1.5])
        line(xlim, [-1.5 -1.5])
        
        line(xlim, [2 2])
        line(xlim, [-2 -2])
        
        % baseline
        WNonday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        line([WNonday-0.5 WNonday-0.5], ylim, 'Color' ,'r');
        
        % day I called consol start
        try
        day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
        day2=day1+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins-1;
        
        line([day1-0.5 day1-0.5], ylim);
        line([day2+0.5 day2+0.5], ylim);
        catch err
        end
    end
end










            
            
    



