function [Params, RawDatStruct]=lt_seq_dep_pitch_DayRawDat(Params, plotON, saveON, phrase, plotLMANinact)
%% 6/15/15 - added catch soing.
%% LT 5/3/15 - added stuff to look at muscimol in LMAN data -
% i.e. plotting raw data over time, variability, contours.

% plotLMANinact = 1; plots. keep as 0 to not plot, if don't define, then will replace with whatever plotON is.

% example:
% Params.DayRawDat.fs=32000;
% Params.DayRawDat.pc_harms=1; % harmonics to take weighted avg over. 1 or 2 is good.
% Params.DayRawDat.batch='batch.labeled.all';
% Params.DayRawDat.syllables={'a','b','c'};
% Params.DayRawDat.frequency_range={[1300 2200], [2800 3950],[2150 3150]};
% Params.DayRawDat.pc_dur=[0.12,0.09,0.11];
% Params.DayRawDat.pc_time_window={[375 525],[60 220],[55 320]};
% Params.DayRawDat.pc_sigma=1;
%
% % plot and save?
% plotON=1;
% saveON=1;
%
% % Related to LMAN inactivation
% plotLMANinact=1;
% Params.DayRawDat.Musc_On_Time='1400'; % time given muscimol - will plot data with temporal lag after this.


% phrase = '' or empty - fills automatically with day phrase (used for
% deciding what folder to save in.



%% LT 4/28/15 - added catch trial as 11th column

%% LT- 1/31/15 modified from lt_compile_RawDatStruct_data

%% LT - 11/13/14 -
% Compiles pitch data - pitch contour, summary pitch data, and sequence
% data. Also
% removes outliers based on the pitch contour (i.e. any timebin with
% extreme values removes the entire syllable data pt).

% Once this is run, you can use lt_compile_RawDatStruct_data_SEQFILTER to
% filter the data in a sequence dependent
% manner (e.g. want not just b, but b that follows "a" specifically.
% Run in day folder.  Makes structure and saves in a folder one dir up.

% examples inputs:
% batch='batch.labeled.all';
% syllables={'b','c'};
% frequency_range={[3000 3900],[2000 3000]}; % for findwnote
% pc_dur=[0.06,0.1]; % size of syl data window (in sec); (how much data to get for each rend (relative to onset), in sec)
% pc_time_window={[25 120],[45 220]}; for pitch contour (time bins to avg).
% SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should combine
% SylTargList={'c','b'}; % these must already have raw data compiled above
% SeqPostList={'b',''};
% plotON=1;


% PARAMETERS
% findwnote
% calcFFT=0; % don't get FFT at arbitrary timeshift
% fs=32000;
%
% % pitch contour
% pc_harms=1; % harmonics to take weighted avg over. 1 or 2 is good.
% plotON=1;
%
% % AUTO params
% NumSyls=length(syllables);
% [birdname bluejaynum date phrase]=lt_get_birdname_date_from_dir(1);
% timestamp=lt_get_timestamp;
% curr_dir=pwd;


%% Params

if ~exist('plotLMANinact','var');
    plotLMANinact=0;
end

if plotLMANinact==1;
    plotON=1;
end

% 1) Take input params and simplify variable name
fs=Params.DayRawDat.fs;
batch=Params.DayRawDat.batch;
syllables=Params.DayRawDat.syllables;

% pitch contour
pc_harms=Params.DayRawDat.pc_harms; % harmonics to take weighted avg over. 1 or 2 is good.
pc_sigma=Params.DayRawDat.pc_sigma; % std of gaussian kernel for spectrogram
frequency_range=Params.DayRawDat.frequency_range;

% dat_dur=Params.DayRawDat.dat_dur; % duration of data to collect
pc_dur=Params.DayRawDat.pc_dur;
pc_time_window=Params.DayRawDat.pc_time_window;


% 2) AUTO params
NumSyls=length(Params.DayRawDat.syllables);
[birdname bluejaynum date phrase_tmp]=lt_get_birdname_date_from_dir(1);
timestamp=lt_get_timestamp;
curr_dir=pwd;

if ~exist('phrase','var') || isempty(phrase);
    phrase = phrase_tmp;
end


%% FIRST, collect syllable data - sound data, pitch contour, spectrogram, timestamps, sequence
SylsToDrop=[]; % will fill up with syls to skip (no data);
NumSyls2=NumSyls;

for i=1:NumSyls;
    try % i.e. maybe this day does not have this syllable
        fvalsstr_forpc_all.(syllables{i}) =...
            findwnote2tw_v4_LT(batch, syllables{i},'',-0.005,...
            frequency_range{i},pc_dur(i)*fs,1,1,'obs0',0);
        
        NN=length(fvalsstr_forpc_all.(syllables{i}));
        
        % PUT INTO OUTPUT STRUCTURE
        for ii=1:NN; % all renditions
            RawDatStruct.data.(syllables{i}){ii,3}=fvalsstr_forpc_all.(syllables{i})(ii).datt;
            RawDatStruct.data.(syllables{i}){ii,5}=fvalsstr_forpc_all.(syllables{i})(ii).fn;
            RawDatStruct.data.(syllables{i}){ii,6}=fvalsstr_forpc_all.(syllables{i})(ii).datenum;
            RawDatStruct.data.(syllables{i}){ii,7}=fvalsstr_forpc_all.(syllables{i})(ii).lbl;
            RawDatStruct.data.(syllables{i}){ii,8}=fvalsstr_forpc_all.(syllables{i})(ii).NotePos;
            RawDatStruct.data.(syllables{i}){ii,9}=fvalsstr_forpc_all.(syllables{i})(ii).TRIG;
            RawDatStruct.data.(syllables{i}){ii,11}=fvalsstr_forpc_all.(syllables{i})(ii).CatchTrial;
            RawDatStruct.data.(syllables{i}){ii,13}=fvalsstr_forpc_all.(syllables{i})(ii).CatchSong;
            
            x=fvalsstr_forpc_all.(syllables{i})(ii).NotePos;
            RawDatStruct.data.(syllables{i}){ii,10}=fvalsstr_forpc_all.(syllables{i})(ii).offs(x)...
                -fvalsstr_forpc_all.(syllables{i})(ii).ons(x);
        end
        
    catch err % if syl does not exist, then remove it from analysis.
        NumSyls2=NumSyls2-1;
        SylsToDrop=[SylsToDrop i];
    end
    
end

syllables(SylsToDrop)=[]; % removing NoData syls.
frequency_range(SylsToDrop)=[];
pc_dur(SylsToDrop)=[];
pc_time_window(SylsToDrop)=[];

NumSyls=NumSyls2;

% how many renditions for each syl.
for i=1:NumSyls;
    NumRends(i)=length(fvalsstr_forpc_all.(syllables{i}));
end

% Update params
Params.DayRawDat.syllables=syllables;
Params.DayRawDat.frequency_range=frequency_range;
Params.DayRawDat.pc_dur=pc_dur;
Params.DayRawDat.pc_time_window=pc_time_window;



%% Second, COLLECT PITCH CONTOUR DATA

for i=1:NumSyls;
    [PC_raw.(syllables{i}), pc_F{i}, pc_T{i}, ~]=jc_pitchcontourFV_LT(fvalsstr_forpc_all.(syllables{i}),...
        1024,1020,1, frequency_range{i}(1),...
        frequency_range{i}(2),pc_harms,'obs0');
    
    for ii=1:NumRends(i); % all renditions
        RawDatStruct.data.(syllables{i}){ii,2}=PC_raw.(syllables{i})(:,ii)'; % pc
        % removed spectrogram becuase takes too much space - for 300 syls,
        % about 30MB
        %         RawDatStruct.data.(syllables{i}){ii,4}=spec{ii}; % spectrogram
        
    end
end


%% EXTRACT FF DATA FROM PITCH CONTOUR

for i=1:length(syllables);
    PCmat=cell2mat(RawDatStruct.data.(syllables{i})(:,2)); % cell to mat
    PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
    
    N=size(RawDatStruct.data.(syllables{i}),1); % num rends
    RawDatStruct.data.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
end



%% Assign legend (what do data columns mean?)
Params.DayRawDat.Legend_data{1,:}={'1_FF','2_PitchContour','3_SoundDat','4_Spec',...
    '5_FileName','6_DateNum','7_LabelStr','8_NotePos','9_Trig','10_NoteDur','11_CatchTrial', '13_CatchSong'};


%% PLOT
if plotON==1;
    % PLOT PITCH CONTOUR
    for i = 1:NumSyls;
        %plots the mean and individual renditions
        figure; hold on;
        subplot(2,1,1), hold on;
        plot(PC_raw.(syllables{i}),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(mean(PC_raw.(syllables{i})'),'Linewidth',2)
        
        % plot lines with bins vals
        line([pc_time_window{i}(1) pc_time_window{i}(1)],ylim);
        line([pc_time_window{i}(2) pc_time_window{i}(2)],ylim);
        
        xlabel(['Time, bins: ' num2str(pc_T{i}(2)-pc_T{i}(1))])
        ylabel('Frequency (Hz)')
        title(['Mean and individual pitch contours for ' syllables{i}]);
        
        % plot with real time axis
        subplot(2,1,2), hold on;
        plot(pc_T{i}*1000,PC_raw.(syllables{i}),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(pc_T{i}*1000,mean(PC_raw.(syllables{i})'),'Linewidth',2)
        xlabel('Time (ms)')
        ylabel('Frequency (Hz)')
        title(['Mean and individual pitch contours for ' syllables{i}]);
        line([pc_T{i}(pc_time_window{i}(1))*1000 pc_T{i}(pc_time_window{i}(1))*1000],ylim)
        line([pc_T{i}(pc_time_window{i}(2))*1000 pc_T{i}(pc_time_window{i}(2))*1000],ylim)
    end
    
    % PLOT FF values
    %     figure; hold on;
    %     title('Individual renditions, FF');
    
    
    
    % for i=1:length(syllables);
    %     PCmat=cell2mat(RawDatStruct.data_WithOutlier.(syllables{i})(:,2)); % cell to mat
    %     PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
    %
    %     N=size(RawDatStruct.data_WithOutlier.(syllables{i}),1); % num rends
    %     RawDatStruct.data_WithOutlier.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
    % end
    %
    %
    
end
%% ---- PLOT THINGS USEFUL FOR LOOKING AT LMAN INACTIVATION
% to sort different conditions (E.g. PBS, MUSC), put those words in the
% song file name (e.g. pu11wh87_MUSC_030515_134138.71.cbin);

RunBin=20; % 10 renditions
plot_cols=lt_make_plot_colors(NumSyls,0,0);

ExptCondition_codes={'PBS','MUSC'}; % i.e. if filename says "PBS", then is code 1 (musc=2 and so forth)

if plotLMANinact==1;
    for i = 1:NumSyls;
        syl=syllables{i};
        
        % ==== SORT TODAY's TRIALS INTO CONDITIONS ==================
        % get all filenames
        filenames=RawDatStruct.data.(syl)(:,5);
        % get conditions from filenames
        ExptConditions=[];
        for ii=1:length(filenames);
            underscores=strfind(filenames{ii},'_');
            exptcond=filenames{ii}(underscores(1)+1:underscores(2)-1);
            
            % convert that condition (string) to code (1, 2, ...);
            switch exptcond
                case 'PBS';
                    ExptConditions(ii)=1;
                case 'MUSC';
                    ExptConditions(ii)=2;
                otherwise
                    disp('PROBLEM - a song is not PBS or MUSC based on filename');
            end
        end
        
        
        % -- SAVE EXPERIMENTAL CONDITION TO OUTPUT STRUCTURE --------
        RawDatStruct.data.(syl)(:,12)=num2cell(ExptConditions');
        
        % update legend in params
        Params.DayRawDat.Legend_data{1,:}={'1_FF','2_PitchContour','3_SoundDat','4_Spec',...
            '5_FileName','6_DateNum','7_LabelStr','8_NotePos','9_Trig','10_NoteDur','11_CatchTrial','12_ExptCond'};
        Params.DayRawDat.ExptCondition_codes=ExptCondition_codes;
        
    end
    
    %
    %     % ==== DOES TODAY HAVE PBS/MUSC, OR BOTH? ==================
    %
    %     if any(num2cell(ExptConditions')==1); % then PBS exists
    %
    %
    %
    
    
    % ==== PLOT VARIOUS THINGS ==============================
    plotcols={[0 0 1],[1 0 0]};
    for i = 1:NumSyls;
        syl=syllables{i};
        
        lt_figure; hold on;
        
        % get pitch and time and expt cond values (sorted into condition)
        FFvals=cell2mat(RawDatStruct.data.(syl)(:,1));
        
        [~, tmp]=lt_convert_datenum_to_hour(cell2mat(RawDatStruct.data.(syl)(:,6)));
        Tvals=tmp.hours;
        
        ExptCondvals=cell2mat(RawDatStruct.data.(syl)(:,12));
        
        PCmat=cell2mat(RawDatStruct.data.(syl)(:,2));
        
        % -- Sort all those values by time
        [Tvals, inds]=sort(Tvals);
        FFvals=FFvals(inds);
        ExptCondvals=ExptCondvals(inds);
        PCmat=PCmat(inds,:);
        
        
        % 1) -- plot FF vs. time (all values)
        lt_subplot(2,2,1); hold on;
        title('FF vs. time'); xlabel('time (hrs)');
        ylabel('FF (hz)');
        
        % Plot
        % PBS
        lt_plot(Tvals(ExptCondvals==1),FFvals(ExptCondvals==1),{'Color','b','MarkerSize',4});
        
        % MUSC
        lt_plot(Tvals(ExptCondvals==2),FFvals(ExptCondvals==2),{'Color','r','MarkerSize',4});
        
        % -- plot lines to indicate when muscimol started and ended
        % (optional)
        if isfield(Params.DayRawDat,'Musc_On_Time');
            % convert from string (e.g. '1432' to hours)
            tmp=datenum(Params.DayRawDat.Musc_On_Time,'HHMM');
            [~, ontime]=lt_convert_datenum_to_hour(tmp);
            
            line([ontime.hours ontime.hours],ylim,'Color','r','LineWidth',2);
        end
        if isfield(Params.DayRawDat,'Musc_Off_Time');
            % convert from string (e.g. '1432' to hours)
            tmp=datenum(Params.DayRawDat.Musc_Off_Time,'HHMM');
            [~, offtime]=lt_convert_datenum_to_hour(tmp);
            
            line([offtime.hours offtime.hours],ylim,'Color','r','LineWidth',2);
        end
        
        
        
        % 2)  -- Plot running FF + STD
        %         subplot(2,2,1); hold on;
        %         % -- plot running avg + running std
        %         title(['Running mean +/- std. Binsize: ' num2str(RunBin)]);
        %         xlabel('time (hrs)');
        %         ylabel('FF (hz)');
        
        % calculate running stats
        FFvals_sm=lt_running_stats(FFvals,RunBin);
        Tvals_sm=lt_running_stats(Tvals,RunBin);
        
        % Plot
        if length(Tvals_sm.Mean)>1;
        shadedErrorBar(Tvals_sm.Mean,FFvals_sm.Mean,FFvals_sm.STD,{'Color','k','LineWidth',2},1);
        end
        
        
        % 2) --- FF vs. rendition number
        lt_subplot(2,2,2); hold on;
        title('FF vs. rendition'); xlabel('rendition num');
        ylabel('FF (hz)');
        
        Xvals=1:length(Tvals);
        
        % Plot
        % PBS
        lt_plot(Xvals(ExptCondvals==1),FFvals(ExptCondvals==1),{'Color','b','MarkerSize',4});
        
        % MUSC
        lt_plot(Xvals(ExptCondvals==2),FFvals(ExptCondvals==2),{'Color','r','MarkerSize',4});
        
        % calculate running stats
        FFvals_sm=lt_running_stats(FFvals,RunBin);
        Xvals_sm=lt_running_stats(Xvals,RunBin);
        
        % Plot
                if length(Xvals_sm.Mean)>1;
        shadedErrorBar(Xvals_sm.Mean,FFvals_sm.Mean,FFvals_sm.STD,{'Color','k','LineWidth',2},1);
                end
        
        
        % 3) -- Running CV vs time
        lt_subplot(2,2,3); hold on;
        title(['Running CV, plotted at mean timepoint; Binsize: ' num2str(RunBin)]);
        xlabel('time (hrs)');
        ylabel('CV');
        
        CV=FFvals_sm.STD./FFvals_sm.Mean;
        
        % plot
        lt_plot(Tvals_sm.Mean,CV);
        
        % -- plot lines to indicate when muscimol started and ended
        % (optional)
        if isfield(Params.DayRawDat,'Musc_On_Time');
            % convert from string (e.g. '1432' to hours)
            tmp=datenum(Params.DayRawDat.Musc_On_Time,'HHMM');
            [~, ontime]=lt_convert_datenum_to_hour(tmp);
            
            line([ontime.hours ontime.hours],ylim,'Color','r','LineWidth',2);
        end
        if isfield(Params.DayRawDat,'Musc_Off_Time');
            % convert from string (e.g. '1432' to hours)
            tmp=datenum(Params.DayRawDat.Musc_Off_Time,'HHMM');
            [~, offtime]=lt_convert_datenum_to_hour(tmp);
            
            line([offtime.hours offtime.hours],ylim,'Color','r','LineWidth',2);
        end
        
        
        
        % 4) -- Running CV vs. rendition
        lt_subplot(2,2,4); hold on;
        title(['Running CV, versus rendition); Binsize: ' num2str(RunBin)]);
        xlabel('Rendition number (mid-window)');
        ylabel('CV');
        
        CV=FFvals_sm.STD./FFvals_sm.Mean;
        
        % plot
        lt_plot(Xvals_sm.Mean,CV);
        
        
        
        
        % 4) -- Wiggle of pitch contour
        % -- calculate deviation from running avg.
        %         PC_run_bins_list=[5 10 15]; % bins to try for running avg, in bins (note around 0.125ms per time bin);
        %
        %         PCmat=cell2mat(RawDatStruct.data.(syl)(:,2)); % d1=trials, d2 = time
        %         PCmat_windowed=PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)); % gets PC for all rends within the desired time window
        %
        %         % within that window, get running avg across time for various
        %         % sized bins
        %         wiggle_integrated_dev=[];
        %         for j=1:size(PCmat,1);
        %             tmp=smooth(PCmat(j,:),20)';
        %
        %             wiggle_integrated_dev(j)=sum(abs(PCmat(j,:)-tmp));
        %
        %         end
        %
        %         subplot(2,2,4); hold on;
        %         title('Wiggle (sum deviation from running avg)');
        %
        %         lt_plot(Tvals,wiggle_integrated_dev);
        
        
        % 5) -- syllable durations
        
        
        % END) label figure
        lt_subtitle(['Syllable: ' syl]);
        
        
        % ==== PLOT MEANS COMPARING MUSC TO PBS EPOCHS
        lt_figure; hold on;
        hsplot=[];
        for ii=1:length(ExptCondition_codes);
            
            vals=FFvals(ExptCondvals==ii); % vals in this condition group
            vals_PC=PCmat(ExptCondvals==ii,:);
            
            if ~isempty(vals); % only continue if have data for this condition.
                
                % if this is MUSC (after PBS) data, and if I have entered start time in params, then
                % throw out the first 45 minutes of data
                notifier=[];
                if strcmp(ExptCondition_codes{ii},'MUSC')
                    if  isfield(Params.DayRawDat, 'Musc_On_Time');
                        % find data that is within 45 minutes from onset
                        
                        % -- ontime in units of hours
                        tmp=datenum(Params.DayRawDat.Musc_On_Time,'HHMM');
                        [~, ontime]=lt_convert_datenum_to_hour(tmp);
                        
                        
                        % Time vals
                        Tvals_MUSC=Tvals(ExptCondvals==ii);
                        
                        % which vals pass criterion?
                        vals=vals(Tvals_MUSC>ontime.hours+45/60);
                        vals_PC=vals_PC(Tvals_MUSC>ontime.hours+45/60,:);
                        
                        % notifier for plotting title
                        notifier=1;
                    end
                end
                
                SUMMARYDATA{ii}.FFmean=mean(vals);
                SUMMARYDATA{ii}.FFstd=std(vals);
                SUMMARYDATA{ii}.N=length(vals);
                SUMMARYDATA{ii}.FFsem=std(vals)/sqrt(length(vals)-1);
                SUMMARYDATA{ii}.PCmean=mean(vals_PC,1);
                SUMMARYDATA{ii}.PCstd=std(vals_PC,0,1);
                SUMMARYDATA{ii}.PCsem=std(vals_PC,0,1)/sqrt(length(vals)-1);
                
                % bootstrap stats for CV
                CV_btstats=lt_bootstrap(vals,'cv',1000);
                SUMMARYDATA{ii}.FFcv_btstrapCI=CV_btstats.CI;
                SUMMARYDATA{ii}.FFcv=CV_btstats.MEAN;
                
                
                % plot mean with sem
                lt_subplot(2,3,1); hold on;
                if notifier==1;
                    title('Mean (SEM) (1st 45 minutes removed)');
                else
                    title('Mean (SEM)');
                end
                
                errorbar(ii,SUMMARYDATA{ii}.FFmean,SUMMARYDATA{ii}.FFsem,...
                    'o','Color','k','MarkerFaceColor','k');
                xlim([-0.5 3.5])
                set(gca,'Xtick',[1 2]);
                set(gca,'XtickLabel',ExptCondition_codes);
                
                % plot CV
                lt_subplot(2,3,2); hold on;
                if notifier==1;
                    title('CV (95%-ile) (1st 45 minutes removed)');
                else
                    title('CV (95%-ile)');
                end
                
                errorbar(ii,SUMMARYDATA{ii}.FFcv, SUMMARYDATA{ii}.FFcv-SUMMARYDATA{ii}.FFcv_btstrapCI(1),...
                    SUMMARYDATA{ii}.FFcv_btstrapCI(2)-SUMMARYDATA{ii}.FFcv,'o','Color','k','MarkerFaceColor','k');
                xlim([-0.5 3.5])
                
                set(gca,'Xtick',[1 2]);
                set(gca,'XtickLabel',ExptCondition_codes);
                
                
                % -- plot distribution of pitch values
                lt_subplot(2,3,3); hold on;
                title('Histogram of FF vals');
                
                [n, centers]=hist(vals,15);
                % normalize to get psd
                n=n./sum(n);
                
                plot(centers,n,'-','Color',plotcols{ii},'LineWidth',2);
                
                ylabel('prob density');
                xlabel('FF (hz)');
                
                
                % -- plot subset of pitch contours + mean
                hsplot(ii)=lt_subplot(2,3,3+ii); hold on;
                title(['random 25 PCs (' ExptCondition_codes{ii}]);
                % plot subset
                % get random subset of PCs
                n=size(vals_PC,1);
                inds=randperm(n);
                
                try
                    inds=inds(1:25);
                catch err
                    disp('Not enough trials to get 25 PCs, using all you got');
                    
                end
                
                plot(vals_PC(inds,:)');
                
                
                % -- overlay mean pitch contours for both conditions over
                % each other.
                lt_subplot(2,3,6); hold on;
                title('Mean (STD) pitch contour');
                
                if length(SUMMARYDATA{ii}.PCmean)>1;
                shadedErrorBar(1:length(SUMMARYDATA{ii}.PCmean),SUMMARYDATA{ii}.PCmean,SUMMARYDATA{ii}.PCstd,{'Color',plotcols{ii},'LineWidth',1.5},1);
                end
            end
        end
        linkaxes(hsplot,'xy');
        
        
        lt_subtitle(['Syllable: ' syl]);
        
    end
end


%% GET SUMMARY STATS ACROSS ALL RENDS

for i=1:NumSyls;
    X=cell2mat(RawDatStruct.data.(syllables{i})(:,1)); % FF vals to matrix
    RawDatStruct.summary_stats.(syllables{i}).meanFF=mean(X);
    RawDatStruct.summary_stats.(syllables{i}).medianFF=median(X);
    RawDatStruct.summary_stats.(syllables{i}).sdFF=std(X);
    RawDatStruct.summary_stats.(syllables{i}).n=length(X);
    RawDatStruct.summary_stats.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(RawDatStruct.data.(syllables{i})(:,2)); % rows: rend, cols: timebin
    RawDatStruct.summary_stats.(syllables{i}).pitchcountour_mean=mean(X,1);
    
    %     % get mean spectrogram
    %     X=[];
    %     for ii=1:size(RawDatStruct.data.(syllables{i})(:,4),1); % num rends
    %         if
    %         X(:,:,ii)=spec.(syllables{i}){ii}; % extract individual specs
    %     end
    %     RawDatStruct.summary_stats.(syllables{i}).spec_mean=mean(X,3);
end

% OUTLIER REMOVED
% for i=1:NumSyls;
%     X=cell2mat(RawDatStruct.data_WithOutlier.(syllables{i})(:,1)); % FF vals to matrix
%     RawDatStruct.summary_stats_WithOutlier.(syllables{i}).meanFF=mean(X);
%     RawDatStruct.summary_stats_WithOutlier.(syllables{i}).medianFF=median(X);
%     RawDatStruct.summary_stats_WithOutlier.(syllables{i}).sdFF=std(X);
%     RawDatStruct.summary_stats_WithOutlier.(syllables{i}).n=length(X);
%     RawDatStruct.summary_stats_WithOutlier.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
%
%     % get mean pitch contour
%     X=cell2mat(RawDatStruct.data_WithOutlier.(syllables{i})(:,2)); % rows: rend, cols: timebin
%     RawDatStruct.summary_stats_WithOutlier.(syllables{i}).pitchcountour_mean=mean(X,1);
%
%     % get mean spectrogram
%     %     X=[];
%     %     for ii=1:size(RawDatStruct.data_NoOutlier.(syllables{i})(:,4),1); % num rends
%     %         X(:,:,ii)=RawDatStruct.data_NoOutlier.(syllables{i}){ii,4}; % extract individual specs
%     %     end
%     %     RawDatStruct.summary_stats_NoOutlier.(syllables{i}).spec_mean=mean(X,3);
%
% end
%


%% assign things to params - for saving purposes
% Compile parameters
% syl specific
Params.DayRawDat.pc_F=pc_F;
Params.DayRawDat.pc_T=pc_T;
Params.DayRawDat.pc_time_window=pc_time_window;


% global
Params.DayRawDat.birdname=birdname;
Params.DayRawDat.bluejaynum=bluejaynum;
Params.DayRawDat.date=date{2};
Params.DayRawDat.phrase=phrase;
Params.DayRawDat.timestamp=timestamp;



%% SAVE DATA

if saveON==1
    
    
    % SAVE
    Params.DayRawDat.savedir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/seq_dep_pitch_' phrase];
    Params.DayRawDat.savename_data=['RawDatStruct' '_' date{2}];
    Params.DayRawDat.savename_params=['Params' '_' date{2}];
    
    try
        cd(Params.DayRawDat.savedir);
    catch err % dir does not exist
        mkdir(Params.DayRawDat.savedir);
        cd(Params.DayRawDat.savedir);
    end
    
    save(Params.DayRawDat.savename_data,'RawDatStruct');
    save(Params.DayRawDat.savename_params,'Params');
    
    
    % save figures
    try cd FIGURES;
    catch err
        mkdir FIGURES
        cd FIGURES
    end
    
    mkdir(date{2});
    cd(date{2});
    
    lt_save_all_figs;
    % ---
    
    
    cd(curr_dir)
    
    disp('DONE!');
    
end

