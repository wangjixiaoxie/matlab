function lt_context_PLOT(AllSongsDataMatrix, Params_input, one_switch_a_day)
%% LT 11/15/15 - this is all the plot code in lt_context_CompileAndPlot.
% Params_alldays.NoteToPlot=1;
% Params_alldays.RunBin=10;
%
% Params_alldays.BoundaryTimes={'15Nov2015-0000', '27Nov2015-0000', '29Nov2015-1248'}; % in format of e.g. 05May2014-1423, these are times of switching in experiment (e.g. turning WN off and on, changing pitch contingency, etc)
% Params_alldays.Edge_Num_Rends = 25; % num rends to call "edges"

% FOR PLOTTING SUMMARY EFFECTS OF PROBES (can just leave undefined if not
% using. only uses if one_switch_a_day is on.)
% Params_alldays.BaselinePhase=1; % then phase 1 is baseline (i.e. home -->
% away is just light, without WN)
% Params_alldays.PhasesToPlotProbesFor=[2 3];  % means that should plot
% home--> awayPROBE for these phases


% Params_alldays.throw_out_if_epoch_diff_days=0; % throws out any transitions that overlap with O/N (potentially 2 per O/N)
% one_switch_a_day=1; % AM: home --> away;  PM: away --> home; analysis will focus on a day as a unit of analysis

%%

global Params_alldays
Params_alldays=Params_input;

%% EXTRACT PARAMS
if ~exist('NoteToPlot','var');
    NoteToPlot=Params_alldays.NoteToPlot;
end

if ~exist('RunBin','var');
    RunBin=Params_alldays.RunBin;
end

if ~exist('one_switch_a_day','var');
    one_switch_a_day=0;
end




%% -- FIRST Sort out just datapoints that are detects of desired note
AllSongsData_toplot=AllSongsDataMatrix(AllSongsDataMatrix(:,4)==NoteToPlot,:);

clear AllSongsDataMatrix


%% Make sure data are in chronological order

[~, inds]=sort(AllSongsData_toplot(:,2)); % list of times

AllSongsData_toplot=AllSongsData_toplot(inds,:);


%% Check if there are boundary times. if so, 1) convert to day indices (rel to first day) and 2) annotate each datapoint to know which experimental phase it was in.

if isfield(Params_alldays,'BoundaryTimes');
    
    % convert days to indices
    eventtimes=lt_convert_EventTimes_to_RelTimes(Params_alldays.firstday, Params_alldays.BoundaryTimes);
    Params_alldays.BoundaryTimes_DayUnits=eventtimes.FinalValue; % units of days (e.g. 1.5 is noon 1st day)
    
    
    
    % == ANNOTATE DATA WITH PHASE INFORMATION
    % e.g. phase 1 = start to boundary time 1; phase 2 = boundary time 2 to
    % bt3 ... and so on.
    
    NumPhases = length(Params_alldays.BoundaryTimes_DayUnits) + 1;
    
    % -- For each phase, determine which datapoints are included in it.  then
    % annotate those datapoints in AllSongsDataMatrix
    for i=1:NumPhases;
        
        % 1) save information of the time window for this phase
        if i==1; % first 1st phase,
            Params_alldays.Phases_DayBounds{i}=[0 Params_alldays.BoundaryTimes_DayUnits(1)];
            
        elseif i==NumPhases; % then it is right-bounded by the last datapoint
            last_datapoint_time=AllSongsData_toplot(end,2);
            Params_alldays.Phases_DayBounds{i}=[Params_alldays.BoundaryTimes_DayUnits(end) last_datapoint_time+0.01]; % add 0.01 so includes last datapoint
            
        else
            Params_alldays.Phases_DayBounds{i}=[Params_alldays.BoundaryTimes_DayUnits(i-1) Params_alldays.BoundaryTimes_DayUnits(i)];
        end
        
        
        % 2) Find data within this time window
        lower_bound=Params_alldays.Phases_DayBounds{i}(1);
        upper_bound=Params_alldays.Phases_DayBounds{i}(2);
        
        inds_in_this_phase=AllSongsData_toplot(:,2)>=lower_bound & AllSongsData_toplot(:,2)<upper_bound;
        
        % -- annotate just those datapoints
        AllSongsData_toplot(inds_in_this_phase, 9)=i;
        
    end
end

% == PLOT graph showing phases
lt_figure; hold on; grid on
title('Experimental phases (i.e. boundaries are changes to conditions');
ylabel('Phase');
xlabel('day');

X=AllSongsData_toplot(:, 2); % days
Y=AllSongsData_toplot(:, 9); % phase num

Params_alldays.Phases_DayBounds
plot(X, Y, '-b');

% put dates of transitions
x=0;
y=1;
string=Params_alldays.firstday;
lt_plot_text(x,y,string);
for i=1:length(Params_alldays.BoundaryTimes_DayUnits);
    
    x=Params_alldays.BoundaryTimes_DayUnits(i);
    y=i;
    
    string=Params_alldays.BoundaryTimes{i};
    
    lt_plot_text(x,y,string);
    
end


%% -- Make legend
Params_alldays.Data_Legend={'1_SongNum', '2_TimeOfSong_days','3_TriggerTime_IncludingPreBuffer','4_NoteNum','5_FF','6_Ampl','7_Hit','8_NoteGrougNum', '9_Experimental_Phase'};


%% [OPTIONAL] FILTER DATA IN AD HOC WAYS

if (0)
    AllSongsData_toplot_orig=AllSongsData_toplot;
    % AllSongsDataMatrix=AllSongsDataMatrix_backupOrig;
    AllSongsData_toplot(AllSongsData_toplot(:,5)<3440,:)=[];
end



%% ========= PLOTS

%% HOW MANY RENDITIONS/SONGS PER DAY?
DaysOfSongs=floor(AllSongsData_toplot(:, 2));
NumDays=DaysOfSongs(end);

for i=1:NumDays
    RendsPerDay(i)=sum(DaysOfSongs==i);
    
    % == HOW MANY SONGS PER DAY
    inds=find(DaysOfSongs==i); % inds of all rends that are in this day
    
    SongsPerDay(i)=length(unique(AllSongsData_toplot(inds,1))); % how many unique songs today?
    
end

lt_figure; hold on; title('num renditions per day');
ylabel('num renditions');

lt_plot(1:NumDays,RendsPerDay);

lt_figure; hold on; title('num songs per day');
ylabel('num songs');

lt_plot(1:NumDays,SongsPerDay);

%% Plot all datapoints

if (0)
    % ========== 1) PLOT (by time)
    lt_figure; hold on;
    title(['FF of all detects of notenum: ' num2str(NoteToPlot) ' (by time)']);
    xlabel('time (hr)');
    ylabel('FF (hz)');
    
    % all rends
    X=AllSongsData_toplot(:,2); % time values
    Y=AllSongsData_toplot(:,5); % freq values
    plot(X, Y,'o','MarkerFaceColor','b');
    
    % overlay online hits
    IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
    X=AllSongsData_toplot(IndsOfTriggeredRends,2); % only those triggered
    Y=AllSongsData_toplot(IndsOfTriggeredRends,5);
    plot(X,Y,'or','MarkerFaceColor','r');
    
    % overlay running avg
    X=AllSongsData_toplot(:,2); % time values
    Y=AllSongsData_toplot(:,5); % freq values
    X_sm=lt_running_stats(X,RunBin);
    Y_sm=lt_running_stats(Y,RunBin);
    
    shadedErrorBar(X_sm.Mean,Y_sm.Mean,Y_sm.SEM,{'-k'},1);
    
    % overlay indication of note group
    X=AllSongsData_toplot(:,2); % time values
    Y=AllSongsData_toplot(:,8); % note group
    
    Ylim=ylim;
    plot(X,Ylim(1)+50+Y.*(Ylim(2)-Ylim(1))/2,'-g');
    
    % overlay experiment time boundaries
    plot_time_boundary_lines
    
    
    % =============== PLOT JUST A SINGLE DAY
    day=3;
    lt_figure; hold on;
    title(['FF of all detects of notenum: ' num2str(NoteToPlot) ' (by time)']);
    xlabel('time (hr)');
    ylabel('FF (hz)');
    
    % all rends
    X=AllSongsData_toplot(:,2); % time values
    Y=AllSongsData_toplot(:,5); % freq values
    days=floor(AllSongsData_toplot(:,2));
    % just that day
    X=X(days==day);
    Y=Y(days==day);
    
    plot(X, Y,'o','MarkerFaceColor','b');
    
    % overlay online hits
    IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
    X=AllSongsData_toplot(IndsOfTriggeredRends,2); % only those triggered
    Y=AllSongsData_toplot(IndsOfTriggeredRends,5);
    days=floor(AllSongsData_toplot(IndsOfTriggeredRends,2));
    
    X=X(days==day);
    Y=Y(days==day);
    
    plot(X,Y,'or','MarkerFaceColor','r');
    
    % overlay indication of note group
    X=AllSongsData_toplot(:,2); % time values
    Y=AllSongsData_toplot(:,8); % note group
    days=floor(AllSongsData_toplot(:,2));
    
    X=X(days==day);
    Y=Y(days==day);
    
    Ylim=ylim;
    plot(X,Ylim(1)+50+Y.*(Ylim(2)-Ylim(1))/2,'-g');
    
    % overlay experiment time boundaries
    plot_time_boundary_lines
    
    
    % ============ 2) PLOT by rendition
    lt_figure; hold on;
    title(['FF of all detects of notenum: ' num2str(NoteToPlot) ' by rendition']);
    xlabel('Rendition number');
    ylabel('FF (hz)');
    
    % all rends
    Y=AllSongsData_toplot(:,5); % freq values
    % plot(Y,'--k');
    plot(Y,'o','MarkerFaceColor','b');
    
    % overlay online hits
    IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
    X=1:size(AllSongsData_toplot,1);
    X=X(IndsOfTriggeredRends);
    Y=AllSongsData_toplot(IndsOfTriggeredRends,5);
    
    plot(X,Y,'or','MarkerFaceColor','r');
    
    % overlay running avg
    X=1:size(AllSongsData_toplot,1);
    Y=AllSongsData_toplot(:,5); % freq values
    
    X_sm=lt_running_stats(X,RunBin);
    Y_sm=lt_running_stats(Y,RunBin);
    
    shadedErrorBar(X_sm.Mean,Y_sm.Mean,Y_sm.SEM,{'-k'},1);
    
    
    
    
    % ================= 3) PLOT just the contours for rendition
    lt_figure; hold on;
    title(['Running mean (sem) of FF, notenum: ' num2str(NoteToPlot) ' (binsize): ' num2str(RunBin)]);
    xlabel('Rendition number');
    ylabel('FF (hz)');
    
    % plot running average
    X=1:size(AllSongsData_toplot,1);
    Y=AllSongsData_toplot(:,5); % freq values
    X_sm=lt_running_stats(X,RunBin);
    Y_sm=lt_running_stats(Y,RunBin);
    
    shadedErrorBar(X_sm.Mean,Y_sm.Mean,Y_sm.SEM,{'-k'},1);
    
    % overlay indication of note group
    Y=AllSongsData_toplot(:,8); % note group
    
    Ylim=ylim;
    plot(X,Ylim(1)+50+Y.*(Ylim(2)-Ylim(1))/2,'-g');
    
    
    
    
    
    % ============== 4) PLOT by song number
    lt_figure; hold on;
    title(['FF of all detects of notenum: ' num2str(NoteToPlot) ' (by song number)']);
    xlabel('Song Number');
    ylabel('FF (hz)');
    
    % all rends
    X=AllSongsData_toplot(:,1); % song number
    Y=AllSongsData_toplot(:,5); % freq values
    plot(X, Y,'o','MarkerFaceColor','b');
    
    % overlay online hits
    IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
    X=AllSongsData_toplot(IndsOfTriggeredRends,1); % only those triggered
    Y=AllSongsData_toplot(IndsOfTriggeredRends,5);
    
    plot(X,Y,'or','MarkerFaceColor','r');
    
    % overlay running avg
    X=AllSongsData_toplot(:,1); % time values
    Y=AllSongsData_toplot(:,5); % freq values
    X_sm=lt_running_stats(X,RunBin);
    Y_sm=lt_running_stats(Y,RunBin);
    
    shadedErrorBar(X_sm.Mean,Y_sm.Mean,Y_sm.SEM,{'-k'},1);
    
end


%% == EXTRACT NOTE GROUPS

% what note groups exist?
NoteGroupList=sort(unique(AllSongsData_toplot(:,8)));
Params_alldays.NoteGroupList=NoteGroupList;

% get note group inds
SORTED_DATA.ByNoteGroup=[];

for i=1:length(NoteGroupList);
    ng=NoteGroupList(i);
    
    SORTED_DATA.ByNoteGroup(i).NGnum=ng;
    SORTED_DATA.ByNoteGroup(i).inds=AllSongsData_toplot(:,8)==ng;
end



%% PLOT AVERAGE CONTOUR FOR ALL NOTE GROUPS - QUICK LEARNING?
plot_IndividuallySmoothed=0; % if 0, then won't do plots of individual contours - not worth it
SmBinRends=RunBin;

hfig=[];
hfig(1)=figure; hold on;
lt_plot_format;
hfig(2)=figure; hold on;
lt_plot_format;

plotcols=lt_make_plot_colors(length(NoteGroupList),0,0);

% === PLOT, first filter out data corresponding to NGs
hspot=[];
hplot2=[];
for i=1:length(NoteGroupList);
    % -- Plot each epoch as a separate contour
    % find start and stop inds for each epoch
    
    figure(hfig(1));
    hspot(i)=lt_subplot(1,length(NoteGroupList),i); hold on;
    title(['Note Group: ' num2str(SORTED_DATA.ByNoteGroup(i).NGnum)]);
    ylabel('FF (hz)'); xlabel('Rendition #');
    
    inds=SORTED_DATA.ByNoteGroup(i).inds;
    inds_padded=[0; inds; 0];    % pad beginning and ending to be able to get if starts or ends with epoch
    
    tmp=diff(inds_padded);
    tmp_offsets=tmp(2:end); % this is for finding offsets
    
    SORTED_DATA.ByNoteGroup(i).epoch_start_inds=find(tmp==1);
    SORTED_DATA.ByNoteGroup(i).epoch_end_inds=find(tmp_offsets==-1);
    
    % make sure same num of start and end inds
    if length(SORTED_DATA.ByNoteGroup(i).epoch_start_inds)~=length(SORTED_DATA.ByNoteGroup(i).epoch_end_inds);
        disp('Problem, number of start and end inds dont match - on keyboard mode');
        keyboard
    end
    
    
    % === 1) Plot individual contours for this NG
    MaxNumSylsInEpoch=max(SORTED_DATA.ByNoteGroup(i).epoch_end_inds-SORTED_DATA.ByNoteGroup(i).epoch_start_inds)+1; % max num syls in one epoch
    Ytot=nan(MaxNumSylsInEpoch, length(SORTED_DATA.ByNoteGroup(i).epoch_start_inds)); % dim1 = syl rends; dim 2 = all epochs
    
    for ii=1:length(SORTED_DATA.ByNoteGroup(i).epoch_start_inds);
        startind=SORTED_DATA.ByNoteGroup(i).epoch_start_inds(ii);
        endind=SORTED_DATA.ByNoteGroup(i).epoch_end_inds(ii);
        
        Y=AllSongsData_toplot(startind:endind,5);
        
        plot(Y,'-','Color',rand*([0.1 0.1 1]));
        
        % collect data
        Ytot(1:(endind-startind+1),ii)=Y;
    end
    
    % ==== 2) overlay mean contour
    if size(Ytot,2)>1
        Ymean=nanmean(Ytot,2);
        Ystd=nanstd(Ytot,0,2);
        Ysem=Ystd./sqrt(sum(~isnan(Ytot),2)-1);
        
        shadedErrorBar(1:length(Ymean),Ymean,Ysem,{'-b','LineWidth',2},1);
        
        NumTrials=size(Ytot,2);
        
        
        % === 3) smooth all trials, then take mean contour
        if plot_IndividuallySmoothed==1;
            Ytot_smtrials=nan(MaxNumSylsInEpoch, length(SORTED_DATA.ByNoteGroup(i).epoch_start_inds));
            for ii=1:NumTrials;
                tmp=lt_running_stats(Ytot(:,ii),SmBinRends);
                
                % put back into a matrix
                Ytot_smtrials(1:length(tmp.Mean),ii)=tmp.Mean';
            end
            
            % get mean of smoothed trials
            Ymean_sm=nanmean(Ytot_smtrials,2);
            Ystd_sm=nanstd(Ytot_smtrials,0,2);
            Ysem_sm=Ystd_sm./real(sqrt(sum(~isnan(Ytot_smtrials),2)-1));
            
            plot(1:length(Ymean_sm), Ymean_sm,'-r','LineWidth',2);
        end
        
        % ==== 4) PLOT ONLY mean contours for all NGs (separate fig)
        if plot_IndividuallySmoothed==1; % then first smooths, then  plots mean
            figure(hfig(2));
            title('first smooth, then take mean');
            xlabel('Rendition #');
            ylabel('FF (hz)');
            
            % only take data that is not nan
            Ymean_sm=Ymean_sm(~isnan(Ysem_sm));
            Ysem_sm=Ysem_sm(~isnan(Ysem_sm));
            
            shadedErrorBar(1:length(Ymean_sm), Ymean_sm,Ysem_sm,{'Color',plotcols{i},'LineWidth',2},1);
            text(length(Ymean_sm),Ymean_sm(end),['NoteGroup: ' num2str(NoteGroupList(i)) '; N=' num2str(NumTrials)],'FontSize',12,'Color',plotcols{i});
        else
            % plots mean without smoothing
            figure(hfig(2));
            title('Mean (SEM) of pitch vs. rend, smoothed');
            xlabel('Rendition #');
            ylabel('FF (hz)');
            
            % only take data that is not nan
            Ymean=Ymean(~isnan(Ysem));
            Ysem=Ysem(~isnan(Ysem));
            
            shadedErrorBar(1:length(Ymean), Ymean,Ysem,{'Color',plotcols{i},'LineWidth',2},1);
            text(length(Ymean),Ymean(end),['NoteGroup: ' num2str(NoteGroupList(i)) '; N=' num2str(NumTrials)],'FontSize',12,'Color',plotcols{i});
            
            
            
        end
        % -- Aside -- display sample size
        %     disp(['Number trials for notegroup ' num2str(NoteGroupList(i)) ': ' num2str(NumTrials)])
        
        
    end
end
linkaxes(hspot,'xy');


% % Troubleshooting - confirm inds correspond to note groups
% % -- PLOT by rendition
% figure; hold on;
% title(['FF of all detects of notenum: ' num2str(NoteToPlot) ]);
%
% % all rends
% Y=AllSongsData_toplot(:,5); % freq values
% plot(Y,'-o','MarkerFaceColor','b');
%
% % overlay online hits
% IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
% X=1:size(AllSongsData_toplot,1);
% X=X(IndsOfTriggeredRends);
% Y=AllSongsData_toplot(IndsOfTriggeredRends,5);
%
% plot(X,Y,'or','MarkerFaceColor','r');
%
% % overlay running avg
% X=1:size(AllSongsData_toplot,1);
% Y=AllSongsData_toplot(:,5); % freq values
%
% X_sm=lt_running_stats(X,RunBin);
% Y_sm=lt_running_stats(Y,RunBin);
%
% shadedErrorBar(X_sm.Mean,Y_sm.Mean,Y_sm.SEM,{'-k'},1);
%
% % overlay NG info
% plot(epoch_start_inds,3500,'xy');
% plot(epoch_end_inds,3500,'sy');
%
%
% xlabel('Rendition number');
% ylabel('FF (hz)');


%% PLOT WITH EACH EPOCH AS A DATAPOINT

if ~exist('RendNumBins', 'var');
    RendNumBins=Params_alldays.Edge_Num_Rends; % number of renditions at starts of epochs to analyze. (for analysis of transitions, for instance)
    % RendNumBins=5; % number of renditions at starts of epochs to analyze. (for analysis of transitions, for instance)
end


% == PLOT NUMBER OF RENDITIONS PER EPOCH
lt_figure; hold on;
hplot=[];
for i=1:length(NoteGroupList);
    lt_subplot(1,length(NoteGroupList),i); hold on;
    xlabel('renditions');
    ylabel('count');
    
    hist(SORTED_DATA.ByNoteGroup(i).epoch_end_inds-SORTED_DATA.ByNoteGroup(i).epoch_start_inds+1);
    title(['Note group: ' num2str(NoteGroupList(i))]);
end

lt_subtitle('Histogram of number of renditions per epoch');
% linkaxes(hplot, 'x');


% == COLLECT STATS (one datapoint for each epoch) + PLOT
for i=1:length(NoteGroupList);
    
    % == For each epoch, get stats
    NumEpochs=length(SORTED_DATA.ByNoteGroup(i).epoch_start_inds);
    
    for ii=1:NumEpochs;
        startind=SORTED_DATA.ByNoteGroup(i).epoch_start_inds(ii);
        endind=SORTED_DATA.ByNoteGroup(i).epoch_end_inds(ii);
        
        % -- Collect raw data for this epoch
        FFvals=AllSongsData_toplot(startind:endind,5); % FF
        Tvals=AllSongsData_toplot(startind:endind,2); % time
        SongNumvals=AllSongsData_toplot(startind:endind,1); % song number
        PhaseNumvals = AllSongsData_toplot(startind:endind,9);  % phase number
        
        
        % == Collect summary stats (ENTIRE EPOCH)
        % raw data
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.RAW.FFvals{ii}=FFvals;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.RAW.Tvals{ii}=Tvals;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.RAW.SongNumVals{ii}=SongNumvals;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.RAW.PhaseNumvals{ii}=PhaseNumvals;
        
        
        % summary stats
        if startind>1; % i.e. a previous epoch exists.
            SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_LastSongInPreviousEpoch(ii)=AllSongsData_toplot(startind-1,2); % doesn't matter what epoch number last epoch was.
        else
            SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_LastSongInPreviousEpoch(ii)=Tvals(1); % make this the start of this epoch instead.
        end
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Start(ii)=Tvals(1);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Median(ii)=median(Tvals);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_End(ii)=Tvals(end);
        
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.meanFF(ii)=mean(FFvals);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.medianFF(ii)=median(FFvals);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.stdFF(ii)=std(FFvals);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.semFF(ii)=lt_sem(FFvals);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.N(ii)=length(FFvals);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.NumSongs(ii)=length(unique(SongNumvals));
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.cvFF(ii)=std(FFvals)/mean(FFvals);
        
        % NOTE ON PHASE: want to get one number for each epoch. but each epoch has multiple trials.
        % it is ok, becasue each epoch will usually have all data be in one phase (as I would
        % not change params in middle of epoch. however sometimes changes
        % params, and restarts epoch, giving two phases.
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.PhaseNum{ii}=unique(PhaseNumvals); % this can be multiple numbers, but usually is one.
        
        
        %         bootstats=lt_bootstrap(FFvals,'cv');
        %         SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.cvFF(ii)=bootstats.MEAN;
        %         SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.CI_of_cvFF{ii}=bootstats.CI;
        
        % slope
        [~,~,~,~,~,RegressStats]=lt_regress(FFvals,1:length(FFvals),0);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.slope(ii)=RegressStats.slope;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.intercept(ii)=RegressStats.intercept;
        
        
        % == Collect summary stats (EARLY PART OF EPOCH)
        if length(FFvals)<RendNumBins;
            FFvals_early=FFvals;
            Tvals_early=Tvals;
        else
            FFvals_early=FFvals(1:RendNumBins);
            Tvals_early=Tvals(1:RendNumBins);
        end
        
        % -- Phasenums
        if length(PhaseNumvals)<RendNumBins;
            PhaseNumvals_early=PhaseNumvals;
        else
            PhaseNumvals_early=PhaseNumvals(1:RendNumBins);
        end
        
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.N(ii)=length(FFvals_early);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.meanFF(ii)=mean(FFvals_early);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.medianFF(ii)=median(FFvals_early);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.stdFF(ii)=std(FFvals_early);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.cvFF(ii)=std(FFvals_early)/mean(FFvals_early);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.Tvals{ii}=Tvals_early;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.FFvals{ii}=FFvals_early;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.PhaseNums{ii}=unique(PhaseNumvals_early);
        
        
        
        
        
        %         bootstats=lt_bootstrap(FFvals_early,'cv');
        %         SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.CI_of_cvFF{ii}=bootstats.CI;
        
        % slope
        [~,~,~,~,~,RegressStats]=lt_regress(FFvals_early,1:length(FFvals_early),0);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.slope(ii)=RegressStats.slope;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.intercept(ii)=RegressStats.intercept;
        
        
        % == Collect summary stats (END PART OF EPOCH)
        if length(FFvals)<RendNumBins;
            FFvals_late=FFvals;
            Tvals_late=Tvals;
        else
            FFvals_late=FFvals(end-RendNumBins+1:end);
            Tvals_late=Tvals(end-RendNumBins+1:end);
        end
        
        % -- Phasenums
        if length(PhaseNumvals)<RendNumBins;
            PhaseNumvals_late=PhaseNumvals;
        else
            PhaseNumvals_late=PhaseNumvals(end-RendNumBins+1:end);
        end
        
        
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.N(ii)=length(FFvals_late);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.meanFF(ii)=mean(FFvals_late);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.medianFF(ii)=median(FFvals_late);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.stdFF(ii)=std(FFvals_late);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.cvFF(ii)=std(FFvals_late)/mean(FFvals_late);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.Tvals{ii}=Tvals_late;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.FFvals{ii}=FFvals_late;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.PhaseNums{ii}=unique(PhaseNumvals_late);
        
        
        %         bootstats=lt_bootstrap(FFvals_late,'cv');
        %         SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.cvFF(ii)=bootstats.MEAN;
        %         SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.CI_of_cvFF{ii}=bootstats.CI;
        
        % slope
        [~,~,~,~,~,RegressStats]=lt_regress(FFvals_late,1:length(FFvals_late),0);
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.slope(ii)=RegressStats.slope;
        SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.intercept(ii)=RegressStats.intercept;
    end
end

pause;
close all;

%% == PLOT STATS, ACROSS TIME (ONE EPOCH AS ONE DATAPOINT)
% Time: median of epoch

StatToPlotList={'meanFF'};
% StatToPlotList={'meanFF','slope','intercept','N','cvFF', 'medianFF'};
beginning_or_end=[]; % if 0, beginning, if 1, end, if [] (or if not defined) then entire epoch.

for i=1:length(StatToPlotList);
    StatToPlot=StatToPlotList{i};
    % StatToPlot='meanFF';
    % StatToPlot='meanFF';
    plot_epochs(StatToPlot, beginning_or_end, SORTED_DATA, NoteGroupList);
end

pause;
close all;

%% GET STATS OF TRANSITIONS BETWEEN EPOCHS

SORTED_DATA.ByNoteGroupTransitions=[];
throw_out_if_epoch_diff_days=Params_alldays.throw_out_if_epoch_diff_days; % WORKS - VERIFIED - ideally, if 1, then ignores 1st epoch of the day (i.e. ignores any transitions in whihc one of the two epochs within them clips the ON transition)
%(note: throwing out potantially will throw out 2 transitions per
%overnight)

for i=1:length(NoteGroupList);
    
    % == GET STATS ON TRANSITIONS FROM ALL OTHER NOTE GROUPS
    NumEpochs=length(SORTED_DATA.ByNoteGroup(i).epoch_start_inds); % how many epochs?
    
    % for each epoch, ask what note group came before it.
    for ii=1:NumEpochs;
        
        % things about this note group and epoch
        GlobalInd_OfThisEpoch=SORTED_DATA.ByNoteGroup(i).epoch_start_inds(ii); % index that starts current epoch
        Stats_OfThisNoteGroup=SORTED_DATA.ByNoteGroup(i);
        Epoch_OfThisNoteGroup=ii;
        
        % === CASES WHEN THROW OUT THIS EPOCH
        % is first epoch of data - no preceding epoch
        if GlobalInd_OfThisEpoch==1;
            continue
        end
        
        % =========================== GET INFO ON PRECEDING EPOCH.
        % what is the note group of the index directly previous to this
        % one?
        NGnum_of_previous_epoch=AllSongsData_toplot(GlobalInd_OfThisEpoch-1,8); % this is actual NGnum (0, 1, ...), (not 1, 2, ...
        ind=find(NoteGroupList==NGnum_of_previous_epoch); % index that corresponds to preceding note group
        
        
        % confirm that that is the correct note group (previous)
        if SORTED_DATA.ByNoteGroup(ind).NGnum~=NGnum_of_previous_epoch;
            disp('PROBLEM - notegroup num of previous epoch is not what is expected');
            keyboard
        end
        
        % -- which epoch was it for the preceding note group?
        % first, make sure that the number of end inds is equal to the
        % number of epochs with data
        if length(SORTED_DATA.ByNoteGroup(ind).Stats_OneDataPtPerEpoch.END.meanFF)~=length(SORTED_DATA.ByNoteGroup(ind).epoch_end_inds);
            disp('PROBLEM, number of datapts for epoch end inds is different from N for Stats (one per epoch)');
            keyboard;
        end
        
        % second, use the end inds to find the right epoch
        EpochNum_PreviousNoteGroup=find(SORTED_DATA.ByNoteGroup(ind).epoch_end_inds==GlobalInd_OfThisEpoch-1); % 1st ind of this epoch minus 1 should equal last ind of last epoch
        if length(EpochNum_PreviousNoteGroup)~=1;
            disp('PROBLEM - can"t find epoch num of rpevious notegroup');
        end
        
        % +++++++++ IN PROGRESS
        % All rends of this epoch and the preceding epoch must be within
        % the same day
        if throw_out_if_epoch_diff_days==1;
            TimeStart_ThisEpoch=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Start(ii);
            TimeEnd_ThisEpoch=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_End(ii);
            
            TimeStart_PrevEpoch=SORTED_DATA.ByNoteGroup(ind).Stats_OneDataPtPerEpoch.TIME_Start(EpochNum_PreviousNoteGroup);
            TimeEnd_PrevEpoch=SORTED_DATA.ByNoteGroup(ind).Stats_OneDataPtPerEpoch.TIME_End(EpochNum_PreviousNoteGroup);
            
            % THROW OUT IF EPOCHS ARE NOT ON SAME DAY
            if floor(TimeStart_PrevEpoch)~=floor(TimeEnd_ThisEpoch);
                disp('Threw out epoch becuase not on same day as previous epoch:');
                disp(['time prev start: ' num2str(TimeStart_PrevEpoch) '; time this end: ' num2str(TimeEnd_ThisEpoch) ' (days)']);
                disp(['from ng: ' num2str(NGnum_of_previous_epoch) ' to ng: ' num2str(NoteGroupList(i))]);
                disp(' ');
                
                continue
            end
            %         day_of_current_epoch=floor(AllSongsData_toplot(GlobalInd_OfThisEpoch,2));
        end
        
        
        
        
        % NOTE TO SELF
        % things about this note group and epoch
        GlobalInd_OfThisEpoch=SORTED_DATA.ByNoteGroup(i).epoch_start_inds(ii); % index that starts current epoch
        Stats_OfThisNoteGroup=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch;
        Epoch_OfThisNoteGroup=ii;
        
        % things about the previous note group/epoch
        NGnum_of_previous_epoch;
        Stats_OfPreviousNoteGroup=SORTED_DATA.ByNoteGroup(ind).Stats_OneDataPtPerEpoch;
        EpochNum_PreviousNoteGroup;
        
        
        % ======= Collect data for the transition
        % --- METADATA
        
        % phase (both epochs should be in same phase) [ENTIRE EPOCH RENDS]
        phase_current=Stats_OfThisNoteGroup.PhaseNum{Epoch_OfThisNoteGroup};
        phase_previous=Stats_OfPreviousNoteGroup.PhaseNum{EpochNum_PreviousNoteGroup};
        % only keep this value if they are the same (and are only one
        % value)
        if length(phase_current)==1 & phase_current==phase_previous;
            Edge_phase_EntireEpochRends=phase_current;
        else
            Edge_phase_EntireEpochRends=[];
            disp('note to self - phase mismatch (not a problem)');
        end
        
        % phase (both epochs should be in same phase) [EDGE RENDS]
        phase_current=Stats_OfThisNoteGroup.BEGINNING.PhaseNums{Epoch_OfThisNoteGroup};
        phase_previous=Stats_OfPreviousNoteGroup.END.PhaseNums{EpochNum_PreviousNoteGroup};
        if length(phase_current)==1 & phase_current==phase_previous;
            Edge_phase_EdgeRends=phase_current;
        else
            Edge_phase_EdgeRends=[];
            disp('note to self - phase mismatch (not a problem)');
        end
        
        
        
        
        % -- DIFFERENCE IN EDGE
        % mean
        Edge_meanFF_current=Stats_OfThisNoteGroup.BEGINNING.meanFF(Epoch_OfThisNoteGroup);
        Edge_meanFF_previous =Stats_OfPreviousNoteGroup.END.meanFF(EpochNum_PreviousNoteGroup);
        Edge_meanFF_diff= Edge_meanFF_current-Edge_meanFF_previous;
        
        
        % raw vals
        Edge_rawFF_current=Stats_OfThisNoteGroup.BEGINNING.FFvals{Epoch_OfThisNoteGroup};
        Edge_rawFF_previous=Stats_OfPreviousNoteGroup.END.FFvals{EpochNum_PreviousNoteGroup};
        
        
        % median
        Edge_medianFF_current=Stats_OfThisNoteGroup.BEGINNING.medianFF(Epoch_OfThisNoteGroup);
        Edge_medianFF_previous=Stats_OfPreviousNoteGroup.END.medianFF(EpochNum_PreviousNoteGroup);
        Edge_medianFF_diff= Edge_medianFF_current-Edge_medianFF_previous;
        
        
        % -- DIFFERENCE IN ENTIRE EPOCH
        % cv
        Y_current=Stats_OfThisNoteGroup.cvFF(Epoch_OfThisNoteGroup);
        Y_previous=Stats_OfPreviousNoteGroup.cvFF(EpochNum_PreviousNoteGroup);
        EntireEpoch_cvFF_diff=Y_current-Y_previous;
        
        % slope
        Y_current=Stats_OfThisNoteGroup.slope(Epoch_OfThisNoteGroup);
        Y_previous=Stats_OfPreviousNoteGroup.slope(EpochNum_PreviousNoteGroup);
        EntireEpoch_slopeFF_diff=Y_current-Y_previous;
        
        % mean ff
        Entire_mean_current=Stats_OfThisNoteGroup.meanFF(Epoch_OfThisNoteGroup);
        Entire_mean_previous=Stats_OfPreviousNoteGroup.meanFF(EpochNum_PreviousNoteGroup);
        EntireEpoch_meanFF_diff=Entire_mean_current-Entire_mean_previous;
        
        % median FF
        Entire_median_current=Stats_OfThisNoteGroup.medianFF(Epoch_OfThisNoteGroup);
        Entire_median_previous=Stats_OfPreviousNoteGroup.medianFF(EpochNum_PreviousNoteGroup);
        EntireEpoch_medianFF_diff=Entire_median_current-Entire_median_previous;
        
        
        
        % ============ PUT INTO STRUCTURE
        if ~isfield(SORTED_DATA.ByNoteGroupTransitions, 'NG_first_to_second');
            SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second={};
        end
        
        
        
        % where to put in structure?
        try
            tmp=length(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EpochInds)+1; % if data already exists
        catch err
            SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EpochInds={};
            tmp=1; % start this data
        end
        
        if ~isempty(Edge_phase_EntireEpochRends);
            SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.PhaseNum_BothEpochs(tmp)=Edge_phase_EntireEpochRends;
        else
            SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.PhaseNum_BothEpochs(tmp)=nan;
        end
        
        if ~isempty(Edge_phase_EdgeRends);
            SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.PhaseNum_BothEpochs_EdgeRends(tmp)=Edge_phase_EdgeRends;
        else
            SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.PhaseNum_BothEpochs_EdgeRends(tmp)=nan;
        end
        
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EpochInds{tmp}=[EpochNum_PreviousNoteGroup Epoch_OfThisNoteGroup]; % [epoch of preceding, epochnum of post]
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.StartTime_SecondEpoch(tmp)=Stats_OfThisNoteGroup.TIME_Start(Epoch_OfThisNoteGroup);
        
        
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_rawFF_current{tmp}=Edge_rawFF_current;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_rawFF_previous{tmp}=Edge_rawFF_previous;
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_meanFF_diff(tmp)=Edge_meanFF_diff;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_meanFF_current(tmp)=Edge_meanFF_current;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_meanFF_previous(tmp)=Edge_meanFF_previous;
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_medianFF_diff(tmp)=Edge_medianFF_diff;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_medianFF_current(tmp)=Edge_medianFF_current;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.Edge_medianFF_previous(tmp)=Edge_medianFF_previous;
        
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_cvFF_diff(tmp)=EntireEpoch_cvFF_diff;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_slope_diff(tmp)=EntireEpoch_slopeFF_diff;
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_meanFF_diff(tmp)=EntireEpoch_meanFF_diff;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_meanFF_current(tmp)=Entire_mean_current;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_meanFF_previous(tmp)=Entire_mean_previous;
        
        
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_medianFF_diff(tmp)=EntireEpoch_medianFF_diff;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_medianFF_current(tmp)=Entire_median_current;
        SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ind, i}.EntireEpoch_medianFF_previous(tmp)=Entire_median_previous;
    end
end



%% PLOT STATS, ONE EPOCH TRANSITION AS A DATAPOINT.

RunBin_epochs=5; % smooth over N epochs

% == PLOT ACROSS TIME
% StatToPlot_List={'Edge_meanFF_diff','EntireEpoch_slope_diff', 'EntireEpoch_cvFF_diff', 'EntireEpoch_meanFF_diff'};
StatToPlot_List={'EntireEpoch_meanFF_diff'};
StatToPlot_List={'Edge_meanFF_diff', 'EntireEpoch_meanFF_diff', 'Edge_medianFF_diff', 'EntireEpoch_medianFF_diff'};


for j=1:length(StatToPlot_List);
    StatToPlot=StatToPlot_List{j};
    
    % == PLOT (across time)
    lt_figure; hold on;
    c=1;
    hplot=[];
    
    for i=1:length(NoteGroupList); % the post note group
        
        for ii=1:length(NoteGroupList); % preceding note group
            
            if isempty(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i});
                c=c+1;
                continue
            end
            
            hplot(c)=lt_subplot(length(NoteGroupList),length(NoteGroupList),c); hold on;
            title(['From NG ' num2str(NoteGroupList(ii)) '; To NG:' num2str(NoteGroupList(i))]);
            
            X=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.StartTime_SecondEpoch;  % times of transition (start of epoch 2)
            Y=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.(StatToPlot);  % all diffs
            
            lt_plot(X,Y);
            %             xlabel('Transition time (start of 2nd epoch)');
            %             ylabel([StatToPlot]);
            lt_plot_zeroline;
            
            
            % plot running averages
            if length(X)>RunBin_epochs;
                X_sm=lt_running_stats(X, RunBin_epochs);
                Y_sm=lt_running_stats(Y, RunBin_epochs);
                shadedErrorBar(X_sm.Mean, Y_sm.Mean, Y_sm.SEM, {}, 1);
            end
            
            % plot time boundaries
            plot_time_boundary_lines
            
            c=c+1;
        end
    end
    lt_subtitle(['Y: ' StatToPlot '; X: time at start of epoch2'])
    linkaxes(hplot, 'xy');
end

StatToPlot='Edge_meanFF';
StatToPlot='EntireEpoch_meanFF';

% == PLOT DOTS (ALL DATA) - one figure for each experimental phase
if one_switch_a_day~=1;
    for kk=1:length(Params_alldays.Phases_DayBounds); % num phases
        phasenum=kk;
        phasebounds_days=Params_alldays.Phases_DayBounds{kk};
        
        lt_figure; hold on;
        c=1;
        hplot=[];
        
        for i=1:length(NoteGroupList); % the post note group
            
            for ii=1:length(NoteGroupList); % preceding note group
                
                if isempty(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i});
                    % then this transition did not occur
                    c=c+1;
                    continue
                end
                
                % ===== PLOT
                hplot(c)=lt_subplot(length(NoteGroupList),length(NoteGroupList),c);
                hold on;
                title(['From NG: ' num2str(NoteGroupList(ii)) '; To NG: ' num2str(NoteGroupList(i))]);
                
                % == SORT OUT DATAPOINTS YOU WANT based on phase (i.e. only take
                % those in this phase)
                Phasenum_vals=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.PhaseNum_BothEpochs;
                
                % find the inds that are this phase
                IndsToKeep=Phasenum_vals==phasenum;
                
                % == GET VALUES TO PLOT (phase sorted)
                % --- EDGE MEAN FF
                X=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.([StatToPlot '_previous']); % previous epoch edge
                Y=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.([StatToPlot '_current']); % current edge
                
                % -- sorted
                X=X(IndsToKeep);
                Y=Y(IndsToKeep);
                
                % === skip if x and y are empty
                if isempty(X)
                    c=c+1;
                    continue
                end
                
                % -- PLOT DOT
                plot([1 2], [X' Y'],'-ob');
                xlim([-0.5 3.5])
                
                % -- plot mean
                if length(X)>1;
                    errorbar([1.3 1.7], mean([X' Y']), ([lt_sem(X') lt_sem(Y')]), 'LineStyle','-','Marker','s', 'Color','r','MarkerFaceColor','r');
                end
                
                % == Get p value
                p=signrank(X, Y);
                text(1.7, max(Y)+100, ['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
                
                
                c=c+1;
            end
            
        end
        lt_subtitle([StatToPlot ': Phase num: ' num2str(phasenum) ' (day: ' num2str(Params_alldays.Phases_DayBounds{kk}(1)) ' to ' ...
            num2str(Params_alldays.Phases_DayBounds{kk}(2))]);
        
        linkaxes(hplot, 'xy');
    end
    
    
    % == PLOT DOTS (JUST MEANS /SEM) - one figure for each experimental phase
    for kk=1:length(Params_alldays.Phases_DayBounds); % num phases
        phasenum=kk;
        phasebounds_days=Params_alldays.Phases_DayBounds{kk};
        
        lt_figure; hold on;
        c=1;
        hplot=[];
        
        for i=1:length(NoteGroupList); % the post note group
            
            for ii=1:length(NoteGroupList); % preceding note group
                
                if isempty(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i});
                    % then this transition did not occur
                    c=c+1;
                    continue
                end
                
                % ===== PLOT
                hplot(c)=lt_subplot(length(NoteGroupList),length(NoteGroupList),c);
                hold on;
                title(['From NG: ' num2str(NoteGroupList(ii)) '; To NG: ' num2str(NoteGroupList(i))]);
                
                % == SORT OUT DATAPOINTS YOU WANT based on phase (i.e. only take
                % those in this phase)
                Phasenum_vals=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.PhaseNum_BothEpochs;
                
                % find the inds that are this phase
                IndsToKeep=Phasenum_vals==phasenum;
                
                % == GET VALUES TO PLOT (phase sorted)
                % --- EDGE MEAN FF
                X=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.([StatToPlot '_previous']); % previous epoch edge
                Y=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.([StatToPlot '_current']); % current edge
                
                % -- sorted
                X=X(IndsToKeep);
                Y=Y(IndsToKeep);
                
                % === skip if x and y are empty
                if isempty(X)
                    c=c+1;
                    continue
                end
                
                % -- plot mean
                if length(X)>1;
                    errorbar([1.3 1.7], mean([X' Y']), ([lt_sem(X') lt_sem(Y')]), 'LineStyle','-','Marker','s', 'Color','r','MarkerFaceColor','r');
                end
                xlim([-0.5 3.5])
                
                % == Get p value
                p=signrank(X, Y);
                text(1.7, max(Y)+100, ['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
                
                
                c=c+1;
            end
            
        end
        lt_subtitle([StatToPlot ' Phase num: ' num2str(phasenum) ' (day: ' num2str(Params_alldays.Phases_DayBounds{kk}(1)) ' to ' ...
            num2str(Params_alldays.Phases_DayBounds{kk}(2))]);
        
        linkaxes(hplot, 'xy');
    end
    
    
    
    % == PLOT HISTOGRAM OF DIFFERENCES - one figure for each experimental phase
    StatToPlot='Edge_meanFF';
    StatToPlot='Edge_medianFF';
    % StatToPlot='EntireEpoch_meanFF';
    
    for kk=1:length(Params_alldays.Phases_DayBounds); % num phases
        phasenum=kk;
        phasebounds_days=Params_alldays.Phases_DayBounds{kk};
        
        lt_figure; hold on;
        c=1;
        hplot=[];
        
        for i=1:length(NoteGroupList); % the post note group
            
            for ii=1:length(NoteGroupList); % preceding note group
                
                if isempty(SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i});
                    % then this transition did not occur
                    c=c+1;
                    continue
                end
                
                % ===== PLOT
                hplot(c)=lt_subplot(length(NoteGroupList),length(NoteGroupList),c);
                hold on;
                title(['From NG: ' num2str(NoteGroupList(ii)) '; To NG: ' num2str(NoteGroupList(i))]);
                
                % == SORT OUT DATAPOINTS YOU WANT based on phase (i.e. only take
                % those in this phase)
                Phasenum_vals=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.PhaseNum_BothEpochs;
                
                % find the inds that are this phase
                IndsToKeep=Phasenum_vals==phasenum;
                
                % == GET VALUES TO PLOT (phase sorted)
                % --- EDGE MEAN FF
                Y=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{ii,i}.([StatToPlot '_diff']); % current edge
                
                % -- sorted
                Y=Y(IndsToKeep);
                
                % === skip if y are empty
                if isempty(Y)
                    c=c+1;
                    continue
                end
                
                % -- PLOT histogram
                hist(Y);
                
                % vert line at 0
                line([0 0], ylim, 'Color','k');
                
                % -- plot mean
                if length(Y)>1;
                    line([mean(Y) mean(Y)], ylim, 'Color', 'r');
                end
                
                c=c+1;
            end
            
        end
        lt_subtitle([StatToPlot ': Phase num: ' num2str(phasenum) ' (day: ' num2str(Params_alldays.Phases_DayBounds{kk}(1)) ' to ' ...
            num2str(Params_alldays.Phases_DayBounds{kk}(2))]);
        
        %     linkaxes(hplot, 'xy');
    end
end

pause; close all;

%% ====== PLOTTING DAY BY DAY, for manual context experiments
if one_switch_a_day==1;
    
    lt_figure; hold on;
    xlabel('day (of switch)');
    ylabel('hz');
    title(['Edge mean (binsize=' num2str(Params_alldays.Edge_Num_Rends) ')'])
    numdays=floor(max(AllSongsData_toplot(:,2)));
    
    
    indtmp=find(strcmp(Params_alldays.ConvertToContext1.NoteGroupNum_codes, 'home'));
    ng_home=Params_alldays.ConvertToContext1.NoteGroupNum_codes{indtmp+1};
    
    indtmp=find(strcmp(Params_alldays.ConvertToContext1.NoteGroupNum_codes, 'away'));
    ng_away=Params_alldays.ConvertToContext1.NoteGroupNum_codes{indtmp+1};
    
    indtmp=find(strcmp(Params_alldays.ConvertToContext1.NoteGroupNum_codes, 'awayProbe'));
    ng_awayProbe=Params_alldays.ConvertToContext1.NoteGroupNum_codes{indtmp+1};
    
    plot_epoch_median_times=0;
    
    % ============ HOME --> AWAY
    ng_Ind1=ng_home;
    ng_Ind2=ng_away;
    plotcol='r';
    lt_context_PLOT_sub1;
    
    % ===== AWAY --> HOME
    ng_Ind1=ng_away;
    ng_Ind2=ng_home;
    plotcol='b';
    lt_context_PLOT_sub1;
    
    % ===== HOME ---> AWAY PROBE
    ng_Ind1=ng_home;
    ng_Ind2=ng_awayProbe;
    plotcol='c';
    lt_context_PLOT_sub1;
    
    % ===== AWAY PROBE ---> AWAY
    ng_Ind1=ng_awayProbe;
    ng_Ind2=ng_away;
    plotcol='m';
    lt_context_PLOT_sub1;
    
    
    lt_plot_annotation(1, {'red: HOME-->AWAY', 'blue: AWAY--> HOME', 'cyan: HOME-->AWAYPROBE', ...
        'magen: AWAYPROBE-->AWAY', 'open face=data not contained in day'}, 'k');
    
    plot_time_boundary_lines
    
    
    
    %% ====== COMBINE ALL TRIALS OF SAME TYPE ACROSS DAYS TO GET MEAN [PLOT EDGES]
    
    
    count=1;
    subplotrows=length(Params_alldays.Phases_DayBounds);
    subplotcols=3;
    SubplotsPerFig=subplotrows*subplotcols;
    fignums_alreadyused=[];
    hfigs=[];
    
    hsplots=[];
    
    for i=1:length(Params_alldays.Phases_DayBounds);
        PhaseNum=i;
        
        lt_figure; hold on;
        
        
        % ============ HOME --> AWAY
        ng_Ind1=ng_home;
        ng_Ind2=ng_away;
        plotcol='r';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; home --> away']);
        
        try
            lt_context_PLOT_sub2(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        catch err
        end
        
        % ============ HOME --> AWAY PROBE
        ng_Ind1=ng_home;
        ng_Ind2=ng_awayProbe;
        plotcol='c';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; home --> awayProbe']);
        
        try
            lt_context_PLOT_sub2(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        catch err
        end
        
        % ============ AWAY --> HOME
        ng_Ind1=ng_away;
        ng_Ind2=ng_home;
        plotcol='b';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; away --> home']);
        
        try
            lt_context_PLOT_sub2(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        catch err
        end
        
        
        % === global formatting
        %         lt_subtitle(['Phase: ' num2str(PhaseNum)]);
        
    end
    
    lt_plot_equalaxes(hsplots,0,1);
    
    
    %% ====== COMBINE ALL TRIALS OF SAME TYPE ACROSS DAYS TO GET MEAN [PLOT DIFFERENCES]
    
    
    count=1;
    subplotrows=length(Params_alldays.Phases_DayBounds);
    subplotcols=5;
    SubplotsPerFig=subplotrows*subplotcols;
    fignums_alreadyused=[];
    hfigs=[];
    
    hsplots=[];
    
    for i=1:length(Params_alldays.Phases_DayBounds);
        PhaseNum=i;
        
        
        % ============ HOME --> AWAY
        ng_Ind1=ng_home;
        ng_Ind2=ng_away;
        plotcol='r';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; home->away']);
        
        lt_context_PLOT_sub3(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        
        % ============ HOME --> AWAY PROBE
        ng_Ind1=ng_home;
        ng_Ind2=ng_awayProbe;
        plotcol='c';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; home->awayProbe']);
        
        lt_context_PLOT_sub3(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        
        % ============ AWAY PROBE --> AWAY
        ng_Ind1=ng_awayProbe;
        ng_Ind2=ng_away;
        plotcol='m';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; awayProbe->away']);
        
        lt_context_PLOT_sub3(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        
        
        % ============ AWAY --> awayprobe
        ng_Ind1=ng_away;
        ng_Ind2=ng_awayProbe;
        plotcol='g';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; away->awayProbe']);
        
        lt_context_PLOT_sub3(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        
        % ============ AWAY --> HOME
        ng_Ind1=ng_away;
        ng_Ind2=ng_home;
        plotcol='b';
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        hsplots=[hsplots hsplot];
        title(['Phase: ' num2str(PhaseNum) '; away->home']);
        
        lt_context_PLOT_sub3(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum);
        
    end
    lt_plot_equalaxes(hsplots,0,1);
    
    
    
    
    
    
    %% ====== COMPARE EFFECTS OF PROBE TRIALS (TO BASLINE AND TO EACH OTHER ACROSS PHASES)
    
    
    lt_figure; hold on;
    xval=1;
    Yvals_all={};
    
    % ===== BASELINE (home --> away)
    
    if isfield(Params_alldays, 'BaselinePhase')
        PhaseNum=Params_alldays.BaselinePhase;
        
        % --- home --> away
        ng_Ind1=ng_home;
        ng_Ind2=ng_away;
        plotcol='r';
        textstr='BASELINE';
        
        Yvals=lt_context_PLOT_sub4(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum, xval);
        lt_plot_text(xval, -2+rand*2, textstr, 'b');
        xval=xval+1;
        
        % --- store Yvals for stats
        Yvals_all=[Yvals_all Yvals];
    end
    
    % ===== ALL OTHER PHASES (home --> awayPROBE)
    if isfield(Params_alldays, 'PhasesToPlotProbesFor')
        for i=1:length(Params_alldays.PhasesToPlotProbesFor);
            PhaseNum=Params_alldays.PhasesToPlotProbesFor(i);
            
            % --- home --> awayPROBE
            ng_Ind1=ng_home;
            ng_Ind2=ng_awayProbe;
            plotcol='g';
            textstr=['PH' num2str(PhaseNum)];
            
            Yvals=lt_context_PLOT_sub4(SORTED_DATA, ng_Ind1, ng_Ind2, plotcol, PhaseNum, xval);
            lt_plot_text(xval, -2+rand*2, textstr, 'b');
            xval=xval+1;
            
            % --- store Yvals for stats
            Yvals_all=[Yvals_all Yvals];
        end
        
    end
    
    
    % ===== STATS
    for j=1:length(Yvals_all);
        
        % 1) is each diff from 0?
        [h, p]= ttest(Yvals_all{j});
        
        if p<0.1
            lt_plot_text(j, max(Yvals_all{j})*1.1, ['p=' num2str(p, '%3.2g')], 'r');
        end
        
        % 2) diff from each other?
        for jj=j+1:length(Yvals_all);
            
            [h, p]= ttest2(Yvals_all{j}, Yvals_all{jj});
            
            if p<0.5;
                disp(['bar ' num2str(j) ' vs ' num2str(jj) '; p=' num2str(p, '%3.2g')]);
            end
            
            
            
        end
        
    end
    
    title('BASELINE: home-->away; PHASES: home --> awayProbe');
    
    
end

%% NORMALIZE PROBE (CS+) BY PROBE (CS-) (i.e. subtract CS- mean)
if isfield(Params_alldays, 'Probe_CSplus');
    Probe_CSplus=Params_alldays.Probe_CSplus; % [from to] (actual NG nums)
    Probe_CSminus=Params_alldays.Probe_CSminus;
    
    % Probe_CSplus=[1 2]; % [from to] (actual NG nums)
    % Probe_CSminus=[1 3];
    
    % === DO MULTIPLE STATS
    StatsField={'Edge_meanFF', 'EntireEpoch_medianFF', 'EntireEpoch_meanFF', 'Edge_medianFF'};
    % StatToPlot='EntireEpoch_medianFF';
    
    for kkk=1:length(StatsField);
        StatToPlot=StatsField{kkk};
        
        
        % == convert to indices (1, 2, ...)
        Probe_CSplus_inds(1)=find(NoteGroupList==Probe_CSplus(1));
        Probe_CSplus_inds(2)=find(NoteGroupList==Probe_CSplus(2));
        
        Probe_CSminus_inds(1)=find(NoteGroupList==Probe_CSminus(1));
        Probe_CSminus_inds(2)=find(NoteGroupList==Probe_CSminus(2));
        
        % % == to all data in CSplus, subtract data of CSminus (mean in that phase);
        % for i=1:length(Params_alldays.Phases_DayBounds);
        % phasenum=i;
        %
        % % == Get Probe (CS-) data for that phase
        % Data_ProbeCSminus.phase(i)=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSminus_inds(1), Probe_CSminus_inds(2)};
        %
        % % 1) sort out to get only
        %
        % end
        
        
        % == PLOT - for each phase, plot histogram of differences, after
        % normalizing for CS-
        hfig1=lt_figure; hold on;
        hfig2=lt_figure; hold on;
        hfig3=lt_figure; hold on;
        hplot1=[];
        hplot2=[];
        
        for i=1:length(Params_alldays.Phases_DayBounds);
            phasenum=i;
            
            % ====== FIRST, get mean value for the probe CS-
            Y_ValsToSubtract=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSminus_inds(1), Probe_CSminus_inds(2)}.([StatToPlot '_diff']);
            
            % 1) get only data from this phase num
            phase_nums=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSminus_inds(1), Probe_CSminus_inds(2)}.PhaseNum_BothEpochs;
            
            % 2) keep only data from this phase num
            Y_ValsToSubtract=Y_ValsToSubtract(phase_nums==phasenum);
            
            % 3) Take mean of those values
            Data_ProbeCSminus.phase(phasenum).([StatToPlot '_diff'])=nanmean(Y_ValsToSubtract);
            
            
            % ====== SECOND, plot histogram of differences
            figure(hfig1);
            h=lt_subplot(ceil(length(Params_alldays.Phases_DayBounds)/2), 2, i); hold on; grid on;
            title(['Phase: ' num2str(phasenum)]);
            ylabel(StatToPlot)
            hplot1=[hplot1 h];
            
            
            % 1) collect data, differences during CS+ probe
            Y_diff=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.([StatToPlot '_diff']);
            
            % 2) get only data from this phase num
            phase_nums=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.PhaseNum_BothEpochs;
            
            Y_diff=Y_diff(phase_nums==phasenum);
            
            % 3) subtract same difference from probe cs- trials
            Value_To_Subtract=Data_ProbeCSminus.phase(phasenum).([StatToPlot '_diff']);
            
            Y_diff_norm=Y_diff-Value_To_Subtract;
            
            % 4) Plot histogram
            if any(~isnan(Y_diff_norm)); % only plot if non-nan data exist
                hist(Y_diff_norm);
                xlabel('Pitch difference (Transition to probe epoch)'); ylabel('counts')
            else
                hplot1(end)=[];
                
            end
            lt_subtitle('CS+ probe normalized to CS- probe');
            
            % 5) Plot mean and sem
            line([mean(Y_diff_norm) mean(Y_diff_norm)], ylim, 'LineWidth', 3, 'Color', 'r')
            
            % 6) Save normalized data
            Data_ProbeCSplus.phase(phasenum).([StatToPlot '_diff']).normalized_to_CSminus=Y_diff_norm;
            Data_ProbeCSplus.phase(phasenum).([StatToPlot '_diff']).not_normalized=Y_diff;
            
            
            
            % ====== THIRD, Dot plot
            figure(hfig2);
            h=lt_subplot(ceil(length(Params_alldays.Phases_DayBounds)/2), 2, i); hold on; grid on;
            title(['Phase: ' num2str(phasenum)]); ylabel('Pitch (hz)');
            xlabel('CS+ Probe (PRE --> POST)');
            
            
            hplot2=[hplot2 h];
            
            % 1) collect data
            X=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.([StatToPlot '_previous']);
            Y=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.([StatToPlot '_current']);
            
            % 2) get data just from this phase num
            phase_nums=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.PhaseNum_BothEpochs;
            
            X=X(phase_nums==phasenum);
            Y=Y(phase_nums==phasenum);
            
            % 3) subtract probe difference from Y value
            Y=Y-Data_ProbeCSminus.phase(phasenum).([StatToPlot '_diff']);
            
            % 4) Plot dot plot
            if any(~isnan(Y)); % only plot if non-nan data exist
                % lt_plot(X,Y);
                lt_plot([1 2], [X' Y'], {'LineStyle', '-'});
                %         mean (sem)
                errorbar([1.3 1.7], mean([X' Y']), [lt_sem(X) lt_sem(Y)], 'LineStyle', '-', 'Marker', 's','MarkerSize',7, 'Color', 'r');
                %         sign rank test
                p=signrank(X, Y);
                text(1.5, max(Y)+2, ['p=' num2str(p)],'FontSize', 12, 'FontWeight', 'bold','Color','b');
                
                xlim([-0.5 3.5]);
            else
                hplot2(end)=[];
                
            end
            
            if i==length(Params_alldays.Phases_DayBounds);
                lt_subtitle('CS+ probe normalized to CS- probe');
            end
            
            % === FOURTH, Dot plot, but not normalized
            figure(hfig3);
            
            h=lt_subplot(ceil(length(Params_alldays.Phases_DayBounds)/2), 2, i); hold on; grid on;
            title(['Phase: ' num2str(phasenum)]);
            ylabel(StatToPlot);
            xlabel('CS+ Probe (PRE --> POST)');
            
            % 1) collect data
            X=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.([StatToPlot '_previous']);
            Y=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.([StatToPlot '_current']);
            
            % 2) get data just from this phase num
            phase_nums=SORTED_DATA.ByNoteGroupTransitions.NG_first_to_second{Probe_CSplus_inds(1), Probe_CSplus_inds(2)}.PhaseNum_BothEpochs;
            
            X=X(phase_nums==phasenum);
            Y=Y(phase_nums==phasenum);
            
            % 4) Plot dot plot
            if any(~isnan(Y)); % only plot if non-nan data exist
                lt_plot([1 2], [X' Y'], {'LineStyle', '-'});
                %         mean (sem)
                errorbar([1.3 1.7], mean([X' Y']), [lt_sem(X) lt_sem(Y)], 'LineStyle', '-', 'Marker', 's','MarkerSize',7, 'Color', 'r');
                %         sign rank test
                p=signrank(X, Y);
                text(1.5, max(Y)+2, ['p=' num2str(p)],'FontSize', 12, 'FontWeight', 'bold','Color','b');
                
                xlim([-0.5 3.5]);
            else
                hplot2(end)=[];
                
            end
            if i==length(Params_alldays.Phases_DayBounds);
                lt_subtitle('CS+ probe NOT normalized to CS- probe');
            end
            
        end
    end
    % linkaxes(hplot1, 'xy');
    
    pause;
    if strcmp(input('type y to close all figs' ,'s'),'y');
        close all;
    end
    
    %% ===== COMPARE PROBE EFFECTS (all the difference scores) ACROSS DIFFERENT PHASES
    % what phases to compare?
    phase1=Params_alldays.PhaseToCompare1; % pu35: 4 and 5; rd23: 3 and 4
    phase2=Params_alldays.PhaseToCompare2;
    
    for i=1:length(StatsField);
        StatToPlot=StatsField{i};
        
        lt_figure; hold on;
        
        % == NORMALIZED - Collect data
        hplot=lt_subplot(1,2,1); hold on;
        title('normalized to CS- probe');
        ylabel(StatToPlot);
        xlabel('Phase number')
        
        % phase 1
        Y1=Data_ProbeCSplus.phase(phase1).([StatToPlot '_diff']).normalized_to_CSminus;
        
        % phase 2
        Y2=Data_ProbeCSplus.phase(phase2).([StatToPlot '_diff']).normalized_to_CSminus;
        
        %         % ---
        %     Y1=Y1(1:10);
        %     Y2=Y2(1:10);
        %     % ----
        
        lt_plot(1, Y1');
        lt_plot(2, Y2');
        
        errorbar([1.3 1.7], [mean(Y1) mean(Y2)], [lt_sem(Y1) lt_sem(Y2)], 'LineStyle','none','Color','r','Marker','s', 'MarkerSize',9);
        
        % stats
        p=ranksum(Y1, Y2);
        % [~, p]=ttest2(Y1, Y2);
        
        text(1.5, max(Y1)+1, ['p=' num2str(p)], 'FontSize', 12, 'FontWeight', 'bold','Color','b');
        
        
        xlim([-0.5 3.5]);
        lt_plot_zeroline
        
        set(hplot,'XTick',[1 2]);
        set(hplot,'XTickLabel',{num2str(phase1), num2str(phase2)});
        
        % == NOT NORMALIZED - Collect data
        hplot2=lt_subplot(1,2,2); hold on;
        title('NOT normalized to CS- probe');
        ylabel(StatToPlot);
        xlabel('Phase number')
        
        % phase 1
        Y1=Data_ProbeCSplus.phase(phase1).([StatToPlot '_diff']).not_normalized;
        
        % phase 2
        Y2=Data_ProbeCSplus.phase(phase2).([StatToPlot '_diff']).not_normalized;
        
        %         % ---
        %     Y1=Y1(1:10);
        %     Y2=Y2(1:10);
        %     % ----
        
        
        lt_plot(1, Y1');
        lt_plot(2, Y2');
        
        % stats
        p=ranksum(Y1, Y2);
        % [~, p]=ttest2(Y1, Y2);
        
        text(1.5, max(Y1)+1, ['p=' num2str(p)], 'FontSize', 12, 'FontWeight', 'bold','Color','b');
        
        errorbar([1.3 1.7], [mean(Y1) mean(Y2)], [lt_sem(Y1) lt_sem(Y2)], 'LineStyle','none','Color','r','Marker','s', 'MarkerSize',9);
        
        xlim([-0.5 3.5]);
        lt_plot_zeroline
        
        set(hplot2,'XTick',[1 2]);
        set(hplot2,'XTickLabel',{num2str(phase1), num2str(phase2)});
        
    end
    
end
%% TTEST COMPARING TRANSITION DATA FOR SPECIFIC TRANSTIONS

Transition_one=[1 2]; % [from to] (in Note group num)
Transition_two=[1 3];

% == Convert to indices
Transition_one_indices(1)=find(NoteGroupList==Transition_one(1));
Transition_one_indices(2)=find(NoteGroupList==Transition_one(2));

Transition_two_indices(1)=find(NoteGroupList==Transition_two(1));
Transition_two_indices(2)=find(NoteGroupList==Transition_two(2));


%%
%% ANOVA ON EDGE DIFFERENCE

%% == PLOT THOSE STATS (mean across all epochs).

if (0); % useless, because currently takes average over all experiment - need to separate by dates
    
    lt_figure; hold on;
    for i=1:length(NoteGroupList);
        
        % -- FF mean
        lt_subplot(3,3,1); hold on; title('mean FF');
        ylabel('Hz');
        Y=mean(SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.meanFF);
        Ysem=lt_sem(SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.meanFF);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
        
        % -- FF std
        lt_subplot(3,3,2); hold on; title('std of FF');
        ylabel('Hz');
        Y=mean(SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.stdFF);
        Ysem=lt_sem(SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.stdFF);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
        
        % -- Renditions/song
        lt_subplot(3,3,3); hold on; title('Renditions per song');
        ylabel('Renditions');
        
        Yrends=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.N;
        Ysongs=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.NumSongs;
        
        Y_RendsPerSong=Yrends./Ysongs;
        
        Y=mean(Y_RendsPerSong);
        Ysem=lt_sem(Y_RendsPerSong);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
        
        % -- Epoch duration
        lt_subplot(3,3,4); hold on; title('Epoch duration (first and last song this epoch)');
        ylabel('Hours')
        
        Ydur=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_End...
            -SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Start; % units of days
        % convert to units of hours
        Ydur=Ydur.*24;
        
        Y=mean(Ydur);
        Ysem=lt_sem(Ydur);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
        
        % -- Epoch duration (2nd way)
        lt_subplot(3,3,5); hold on; title('Epoch duration (end last epoch to end this epoch');
        ylabel('Hours')
        
        Ydur=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_End...
            -SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_LastSongInPreviousEpoch; % units of days
        % convert to units of hours
        Ydur=Ydur.*24;
        
        Y=mean(Ydur);
        Ysem=lt_sem(Ydur);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
        
        % -- Epoch duration (2nd way)
        lt_subplot(3,3,5); hold on; title('Epoch duration (end last epoch to end this epoch');
        ylabel('Hours')
        
        Ydur=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_End...
            -SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_LastSongInPreviousEpoch; % units of days
        % convert to units of hours
        Ydur=Ydur.*24;
        
        Y=mean(Ydur);
        Ysem=lt_sem(Ydur);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
        
        % -- Lag from end of last epoch to start of this epoch
        lt_subplot(3,3,6); hold on; title('Time between last epoch (end) and this epoch (start)');
        ylabel('Hours')
        
        Ydur=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Start...
            -SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_LastSongInPreviousEpoch; % units of days
        % convert to units of hours
        Ydur=Ydur.*24;
        
        Y=mean(Ydur);
        Ysem=lt_sem(Ydur);
        errorbar(i, Y, Ysem,'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i},'MarkerSize',7);
        xlim([-0.5, length(NoteGroupList)+1.5])
    end
    
    lt_subtitle('Averages with each epoch contributing one datapoint (+/- SEM using epoch numbers)');
    
end

% == PLOT SAME AS ABOVE, BUT A HISTOGRAM FOR EACH STATISTIC (ACROSS EPOCHS).


% === PLOT ACROSS TIME, each epoch giving one datapoint










%% For a given trial, does its hit status predict the next rendition?

% Plot scattergram - y = next rendition, x = current rendition
% separate data by hits and not hits

HitInds=AllSongsData_toplot(:,7)==1;
NotHitInds=AllSongsData_toplot(:,7)==0;

figure; hold on;
title('Red=hits, Blue=escapes');
ylabel('trials n+1');
xlabel('trials n');

% -- Plot for hits
inds=HitInds;

if any(inds);
    X = AllSongsData_toplot(inds(1:end-1),5);  % trials n (hits)
    X = X-mean(X); % convert to mean subtracted data
    Y = AllSongsData_toplot(logical([0; inds(1:end-1)]),5); % trials n+1
    Y = Y-mean(Y);
    plot(X,Y,'or');
    
    % perform linear regression
    y=Y;
    x=[ones(length(X),1), X];
    
    [b,bint,~,~,stats]=regress(y,x);
    plot(xlim,b(1) + b(2).*xlim,'-r','LineWidth',2);
    plot(xlim,bint(1,1) + bint(2,1).*xlim,'-r');
    plot(xlim,bint(1,2) + bint(2,2).*xlim,'-r');
end



% -- Plot for escapes
inds=NotHitInds;

if any(inds);
    
    X = AllSongsData_toplot(inds(1:end-1),5);  % trials n (hits)
    X = X-mean(X); % convert to mean subtracted data
    Y = AllSongsData_toplot(logical([0; inds(1:end-1)]),5); % trials n+1
    Y = Y-mean(Y);
    plot(X,Y,'ob');
    
    % perform linear regression
    y=Y;
    x=[ones(length(X),1), X];
    
    [b,bint,~,~,stats]=regress(y,x);
    
    plot(xlim,b(1) + b(2).*xlim,'-b','LineWidth',2);
    plot(xlim,bint(1,1) + bint(2,1).*xlim,'-b');
    plot(xlim,bint(1,2) + bint(2,2).*xlim,'-b');
end




% -- This trial hit (or miss) predict how much next trial changes?
figure; hold on;
subplot(1,2,1); hold on;
title('Deviation on next trial (n+1) vs. current trial (n); Red=hits; Blue=escapes');
xlabel('Pitch on trial n (hz)');
ylabel('Deviation (trial n+1 minus n)');


% -- Hit Trials
inds=HitInds;

X = AllSongsData_toplot(inds(1:end-1),5);  % trials n (hits)
X_NextTr = AllSongsData_toplot(logical([0; inds(1:end-1)]),5); % trials n+1

NextTrialDeviation=X_NextTr-X;

plot(X,NextTrialDeviation,'.r');

% perform linear regression
y=NextTrialDeviation;
x=[ones(length(X),1), X];

if any(inds)
    [b,bint,~,~,stats]=regress(y,x);
    plot(xlim,b(1) + b(2).*xlim,'-r','LineWidth',2);
    plot(xlim,bint(1,1) + bint(2,1).*xlim,'-r');
    plot(xlim,bint(1,2) + bint(2,2).*xlim,'-r');
    
    
    % Plot histogram of deviations
    subplot(1,2,2); hold on;
    
    [n, bin]=hist(NextTrialDeviation,30);
    bar(bin,n,'FaceColor',[1 0 0]);
end




% Escapes
subplot(1,2,1); hold on;
inds=NotHitInds;

X = AllSongsData_toplot(inds(1:end-1),5);  % trials n (hits)
X_NextTr = AllSongsData_toplot(logical([0; inds(1:end-1)]),5); % trials n+1

NextTrialDeviation=X_NextTr-X;
plot(X,NextTrialDeviation,'.b');

% perform linear regression
y=NextTrialDeviation;
x=[ones(length(X),1), X];

[b,bint,~,~,stats]=regress(y,x);

plot(xlim,b(1) + b(2).*xlim,'-b','LineWidth',2);
plot(xlim,bint(1,1) + bint(2,1).*xlim,'-b');
plot(xlim,bint(1,2) + bint(2,2).*xlim,'-b');


% Plot histogram of deviations
subplot(1,2,2); hold on;

[n, bin]=hist(NextTrialDeviation,30);
bar(bin,n,'FaceColor',[0 0 1]);


end

%% DOES THE STATISTICS OF HITS WITHIN AN EPOCH EXPLAIN THE TRAJECTORY WITHIN THE EPOCH?
% best looked at with random WN --> tend to hit more lower lead to greater
% positive learning, etc?  Hard to look at with pitch contingent WN,
% because of confound of direction of causality.



%% IS THERE LEARNING AT QUICK TIMESCALE (FEW SONGS).




%% subfunction

function plot_epochs(StatToPlot, beginning_or_end, SORTED_DATA, NoteGroupList)

global Params_alldays

plotcols=lt_make_plot_colors(length(NoteGroupList),0,0);

lt_figure; hold on;
xlabel('Median time of epoch');

hplot=[];
for i=1:length(NoteGroupList);
    
    hplot(i)=lt_subplot(ceil(length(NoteGroupList)/2),2,i); hold on;
    title(['Group number: ' num2str(NoteGroupList(i))]);
    
    % == collect stats to plot
    X=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Median; % time
    
    try
        if beginning_or_end==0; % then beginning data
            Y=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.(StatToPlot);
        elseif beginning_or_end==1; % then end of epoch data
            Y=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.(StatToPlot);
        else % then mean of entire epoch
            Y=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.(StatToPlot);
        end
    catch err % becuase that stat is not defined (for beginning or end)
        disp('PROBLEM - stat does not exist');
    end
    
    % == PLOT
    lt_plot(X, Y);
    
    % plot lines demarcating expt boundaries
    plot_time_boundary_lines;
    
    if any(strcmp(StatToPlot,'slope'));
        lt_plot_zeroline;
    end
    
    
    
    
end
lt_subtitle([StatToPlot ' (each epoch one datapoint)']);

linkaxes(hplot, 'xy');


% === PLOTTING ALL NGs IN ONE PLOT
X={};
Y={};
figure; hold on;
for i=1:length(NoteGroupList);
    NumEpochs=length(SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.(StatToPlot));
    
    % 1) Collect stats across all epochs
    % == prepare cell arrays (one cell for each group num)
    Y{i}=[]; % stat
    X{i}=[]; % median time
    
    % == Collect stats
    % mean time of that epoch
    X{i}=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.TIME_Median;
    
    % data
    try
        if beginning_or_end==0; % then beginning data
            Y{i}=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.BEGINNING.(StatToPlot);
        elseif beginning_or_end==1; % then end of epoch data
            Y{i}=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.END.(StatToPlot);
        else % then mean of entire epoch
            Y{i}=SORTED_DATA.ByNoteGroup(i).Stats_OneDataPtPerEpoch.(StatToPlot);
        end
    catch err % becuase that stat is not defined (for beginning or end)
    end
    
    % 2) plot
    subplot(1,3,1:2); hold on;
    
    % each epoch separate color
    hplot(i)=plot(X{i},Y{i},'o','MarkerFaceColor', plotcols{i}, 'Color', 'k', 'MarkerSize',7);
    
    if i==length(NoteGroupList); % if this is the last NG
        %     % 3) Line connecting all epochs
        %         Xmat=[];
        %         Ymat=[];
        %         for j=1:length(X);
        %             Xmat=[Xmat X{j}];
        %             Ymat=[Ymat Y{j}];
        %         end
        %         % sort by time
        %         [Xmat, inds]=sort(Xmat);
        %         Ymat=Ymat(inds);
        %         plot(Xmat,Ymat,'-k')
        
        legend(hplot,num2str(NoteGroupList));
        
        if any(strcmp(StatToPlot,'slope'));
            lt_plot_zeroline;
        end
    end
    
    % plot lines demarcating expt boundaries
    plot_time_boundary_lines;
    
    % 4) SUMMARY PLOT
    subplot(2,3,3); hold on;
    
    lt_plot(i, Y{i},{'Color',plotcols{i}});
    errorbar(i+0.2, mean(Y{i}), lt_sem(Y{i}),'s','Color',plotcols{i});
    xlim([0 max(NoteGroupList)+1]);
    
    if any(strcmp(StatToPlot,'slope'));
        lt_plot_zeroline;
    end
    
    lt_subtitle([StatToPlot ' (entire epoch)']);
end
end


function plot_time_boundary_lines
global Params_alldays
if isfield(Params_alldays,'BoundaryTimes_DayUnits');
    % plots in units of days
    % overlay experiment time boundaries
    for i=1:length(Params_alldays.BoundaryTimes_DayUnits);
        line([Params_alldays.BoundaryTimes_DayUnits(i) Params_alldays.BoundaryTimes_DayUnits(i)], ylim,'Color','k');
    end
end
end

%% OLD STUFF PUT HERE



% --- IGNORE - OLD WAY FINDING OUT ONLINE TRIGGER FF. NOW I JUST MARK
% DOWN OFFLINE DATAPOINTS AS EITHER REAL TRIGGER (1) OR NOT (0)
% -- Collect online data
% if (0)
%     if ~isempty(AllData(i).data_OnlineTrigs.ttimes);
%         n=size(AllData(i).data_OnlineTrigs.ttimes,1);
%         % figure out FF of online hits. since was not recorded by evtaf, will
%         % get by comparing to offline evtafsim. compare ttimes, and ask which
%         % ttime is closest in offline and online data.
%
%         ff_tmp=[];
%         for j=1:n; % for all datapouints
%             ttime_online=AllData(i).data_OnlineTrigs.ttimes(j);
%             ttime_list_offline=AllData(i).data_OfflineAllDetects(:,1);
%
%             % which offline ttime is closest to online?
%             %             [minval(c), minind]=min(abs(ttime_list_offline-ttime_online));
%             [minval, minind]=min(abs(ttime_list_offline-ttime_online));
%
%             % Note: empirically seems as if offline ttimes are never more
%             % than 5ms diff from online trig times. Given that multiple
%             % trigs never really happen within trig refract of <10ms, will
%             % use 10ms as the value below which I say offline and online
%             % ttimes are the same.  Above this value there are sometimes
%             % very different values, indicating that there was an online
%             % hit that was not detected online.  I throw these out. Low
%             % fraction (2/289; 7/444, on 2 days for rd23gr89). Seems to be related to WN beign on that
%             % sylalble?
%
%             % Troubleshooting - looking at the songs where offline failed
%             % to get all ttimes
%             %             if minval(c)>10;
%             %                 AllData(i).filename
%             %                 AllData(i).data_OnlineTrigs.ttimes
%             %                 j
%             %                 disp(' ');
%             %             end
%
%             if minval<10;
%                 ff_tmp(j)=AllData(i).data_OfflineAllDetects(minind,3); % get FF of the ttime closest
%             else
%                 ff_tmp(j)=nan; % put nan if can't find FF.
%             end
%             %             c=c+1;
%         end
%
%         AllSongsDataMatrix_onlinehits=[AllSongsDataMatrix_onlinehits; hourofsong*ones(n,1), AllData(i).data_OnlineTrigs.ttimes, AllData(i).data_OnlineTrigs.trignotes' ff_tmp'];
%     end
% end
