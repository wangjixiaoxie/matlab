%% LT 3/9/15 - compare two conditions (e.g. stim vs. stim catch) that both have stim timing, looking at effects triggered by stim.

clear OUTSTRUCT

%% NEW PARAMS

Params.STIMTRIG.predur = 60; % how much time (ms) to show precding stim start
Params.STIMTRIG.postdur = 300; % how much time (ms) to show post start of stim



%% AUTO PARAMS
% conversion between PC samples and time (ms)
ms_per_sample_PC=median(diff(1000*Params.tf_bins.tPC));

Params.STIMTRIG.predur_samples=round(Params.STIMTRIG.predur/ms_per_sample_PC);
Params.STIMTRIG.postdur_samples=round(Params.STIMTRIG.postdur/ms_per_sample_PC);

% Extract old params
NumFields=length(Params.FieldsToCheck);
PC_freqrange=Params.PC_freqrange;
pc_harms=Params.pc_harms;
Fs=Params.Fs;
PreDur=1000*Params.PreDur; % time of alignment (i.e. duration of data taken before target syl), ms.
StimDur=Params.StimDur;



%% for each PC, re-zero so that all trials are aligned by stim time

for j=1:NumFields;
    trialfield=Params.FieldsToCheck{j};
    
    NumTrials=size(StatsStruct.(trialfield).PC,2); % d2= trials
    
    for i=1:NumTrials
        
        % Gets stim time for that trial
        stimtime=StatsStruct.(trialfield).TimeSinceLastTrig(i); % stim time in ms preceding t=0;
        stimtime=PreDur-stimtime; % stim time from start of data, in ms
        
        % converts stim time to start and end (in samples)
        [~,stimind]=min(abs(stimtime-1000*Params.tf_bins.tPC)); % finds index corresponding to stimtime
        ind1= stimind-Params.STIMTRIG.predur_samples; % when to start taking data
        ind2= stimind+Params.STIMTRIG.postdur_samples; % when to stop taking data
        
        
        
        % Gets PC, aligned to stim time
        pitchcontour=StatsStruct.(trialfield).PC(:,i); % gets PC for that trial
        
        % check to make sure inds are not negative (start) or past end of
        % data (end)
        if ind1<0;
            disp([trialfield ', trial #: ' num2str(i) ' has negative start ind -skipping']);
            continue
        end
        
        if ind2>length(pitchcontour);
            disp([trialfield ', trial #: ' num2str(i) ' has trial end longer than available data - skipping']);
        continue
        end

        pitchcontour=pitchcontour(ind1:ind2);
        
        OUTSTRUCT.(trialfield).PC_AlignedToStim(:,i)=pitchcontour;
    end
    
end


% Create new t-bins
N=Params.STIMTRIG.predur_samples+Params.STIMTRIG.postdur_samples+1; % number of sampels in each PC
Params.STIMTRIG.tPC=linspace(-Params.STIMTRIG.predur/1000,Params.STIMTRIG.postdur/1000,N);




%% PLOT

figure; hold on;
title('PC aligned to stim time');
PlotCols=lt_make_plot_colors(NumFields,0);


    for i=1:NumFields;
        trialfield=Params.FieldsToCheck{i};
        PC=OUTSTRUCT.(trialfield).PC_AlignedToStim;
        
        NumSyls=size(OUTSTRUCT.(trialfield).PC_AlignedToStim,2);
        
        % Plot all PCs
        subplot(3,1,1); hold on;
        
        plot(1000*Params.STIMTRIG.tPC,PC,'LineStyle','-','Color',PlotCols{i},'LineWidth',0.1) % plot all pitch contours in light shade
        ylabel('Frequency (Hz)')
        title('Individual Pitch Contours');
        xlabel(['Time (ms); aligned at stim times']);
        
        xlim([Params.STIMTRIG.tPC(1)*1000 Params.STIMTRIG.tPC(end)*1000])

        
        % PLOT MEAN AND MEDIAN PC
        subplot(3,1,2); hold on;
        PC_sd=std(PC');
        PC_sem=PC_sd/sqrt(NumSyls-1);
        
        shadedErrorBar(Params.STIMTRIG.tPC*1000,mean(PC'),PC_sem,{'Linewidth',2,'Color',PlotCols{i}},1);
        
        %         % annotate SEM
        %         hfig2(i)=plot(tPC*1000,mean(PC')+PC_sem,'Linewidth',0.5,'Color',PlotCols{i});
        %         plot(tPC*1000,mean(PC')-PC_sem,'Linewidth',0.5,'Color',PlotCols{i});
        
        
        xlim([Params.STIMTRIG.tPC(1)*1000 Params.STIMTRIG.tPC(end)*1000])
        title('Mean +/- SEM  Pitch Contours');
        ylabel('Frequency (Hz)')
        xlim([Params.STIMTRIG.tPC(1)*1000 Params.STIMTRIG.tPC(end)*1000])

        
        % PLOT MEDIAN
        subplot(3,1,3); hold on;
        shadedErrorBar(Params.STIMTRIG.tPC*1000,median(PC'),PC_sem,{'Linewidth',2,'Color',PlotCols{i}},1);
        
        %         % annotate SEM
        %         hfig2(i)=plot(tPC*1000,mean(PC')+PC_sem,'Linewidth',0.5,'Color',PlotCols{i});
        %         plot(tPC*1000,mean(PC')-PC_sem,'Linewidth',0.5,'Color',PlotCols{i});
        
        
        xlim([Params.STIMTRIG.tPC(1)*1000 Params.STIMTRIG.tPC(end)*1000])
        title('Median +/- SEM  Pitch Contours');
        ylabel('Frequency (Hz)')
        xlim([Params.STIMTRIG.tPC(1)*1000 Params.STIMTRIG.tPC(end)*1000])
       
        
        
    end



%% AMPLITUDE PROFILE:

% for amplitude contour, re-zero so that all trials are aligned by stim time

for j=1:NumFields;
    trialfield=Params.FieldsToCheck{j};
    
    NumTrials=size(StatsStruct.(trialfield).PC,2); % d2= trials
    
    for i=1:NumTrials
        
        % Gets stim time for that trial
        stimtime=StatsStruct.(trialfield).TimeSinceLastTrig(i); % stim time in ms preceding t=0;
        stimtime=PreDur-stimtime; % stim time from start of data, in ms
        
        % converts stim time to start and end (in samples)
        [~,stimind]=min(abs(stimtime-1000*Params.tf_bins.tPC)); % finds index corresponding to stimtime
        ind1= stimind-Params.STIMTRIG.predur_samples; % when to start taking data
        ind2= stimind+Params.STIMTRIG.postdur_samples; % when to stop taking data
        
        
        
        % Gets PC, aligned to stim time
        pitchcontour=StatsStruct.(trialfield).PC(:,i); % gets PC for that trial
        
        % check to make sure inds are not negative (start) or past end of
        % data (end)
        if ind1<0;
            disp([trialfield ', trial #: ' num2str(i) ' has negative start ind -skipping']);
            continue
        end
        
        if ind2>length(pitchcontour);
            disp([trialfield ', trial #: ' num2str(i) ' has trial end longer than available data - skipping']);
        continue
        end

        pitchcontour=pitchcontour(ind1:ind2);
        
        OUTSTRUCT.(trialfield).PC_AlignedToStim(:,i)=pitchcontour;
    end
    
end


% Create new t-bins
N=Params.STIMTRIG.predur_samples+Params.STIMTRIG.postdur_samples+1; % number of sampels in each PC
Params.STIMTRIG.tPC=linspace(-Params.STIMTRIG.predur/1000,Params.STIMTRIG.postdur/1000,N);


%% DO SAME FOR AMPLITUDE





