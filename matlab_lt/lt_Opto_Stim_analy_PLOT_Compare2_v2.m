function [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_Compare2(DatStructCompiled,Params)
%% LT 4/30/15 - stopped saving spectrograms, takes up too much memory.

%% LT 4/7/15 - v2, saves overwriting old stuff in directory, instead of making new subdirectory


%% LT 1/12/15 - called by lt_Opto_Stim_analy script
% Given one input structure, and two fields of that structure, will plot
% and compare opto stim song data

% INPUTS:
% DatStructCompiled = DatStructCompiled; % e.g. has fields "All", "Stim", etc, which are sorted
% Params.FieldsToCheck{1}='StimCatch'; % these are field names of that structure. choose two to compare.
% Params.FieldsToCheck{2}='StimNoCatch';

% OUTPUTS:
% StatsStruct - contains processed data.
% Params
% DatStructCompiled

disp('Running...');

%% AUTO PARAMS
NumFields=length(Params.FieldsToCheck);
PC_freqrange=Params.PC_freqrange;
pc_harms=Params.pc_harms;
Fs=Params.Fs;
PreDur=1000*Params.PreDur; % time of alignment (i.e. duration of data taken before target syl), ms.
StimDur=Params.StimDur;

%% Extract Stats for each field

StatsStruct=struct;
fieldstoremove = [];
for k=1:NumFields;
    
    % INITIATE
    fieldname=Params.FieldsToCheck{k};
    
    % make sure has this field
    if isempty(DatStructCompiled.(fieldname))
        fieldstoremove = [fieldstoremove k];
        Params.FieldsToCheck{k}=[];
        continue
    end
    
    %
    DatStruct=DatStructCompiled.(fieldname); % get data structure
    
    NumSyls=length(DatStruct); % how many renditions?
    
    
    % 1) Pitch contour
    
    [StatsStruct.(fieldname).PC, fPC, tPC, ~]=jc_pitchcontourFV_LT(DatStruct,...
        1024,1020,1.5, PC_freqrange(1),...
        PC_freqrange(2),pc_harms,'obs0');
    
    Params.tf_bins.fPC=fPC;
    Params.tf_bins.tPC=tPC;
    
    
    
    % Spectrogram
    for i=1:NumSyls;
        dat=DatStruct(i).datt;
        
        [sm,sp,Params.tf_bins.tSP,Params.tf_bins.fSP]=SmoothData(dat,Fs,1,'‘hanningfirff’'); % ends up doing buttor, with filtfilt.
        
        % REMOVE unwanted frequencies and times
        sp=sp(10:140,:);
        Params.tf_bins.fSP=Params.tf_bins.fSP(10:140);

        % Output Structures
        StatsStruct.(fieldname).sm(:,i)=sm;
        StatsStruct.(fieldname).sp(:,:,i)=sp;

    
        StatsStruct.(fieldname).sm_log(:,i)=log(sm);

    end
    
    
    
    
    % MEAN
    sm_log_mean=mean(StatsStruct.(fieldname).sm_log,2);
    sp_mean=mean(StatsStruct.(fieldname).sp,3);
    
    StatsStruct.(fieldname).sm_log_mean=sm_log_mean;
    StatsStruct.(fieldname).sp_mean=sp_mean;
    
    
    % Entropy - Calculate in sliding window
    sp=StatsStruct.(fieldname).sp;
    for j = 1:size(sp,3); % num renditions
        StatsStruct.(fieldname).WEntropyTimecourse(:,j) = log(geomean(sp(:,:,j))./mean(sp(:,:,j))); % range -inf to 0, verified to be correct LT. for all timepoints.
    end
    
% Time of Day
for     i=1:NumSyls;
    StatsStruct.(fieldname).datenum(i)=DatStructCompiled.(fieldname)(i).datenum;
end


end


% ================== REMOVE EMPTY FIELDS
Params.FieldsToCheck(fieldstoremove) = [];
NumFields = length(Params.FieldsToCheck);


% ASIDE: Adjust timebins so that spectrograms and pitchcountour are aligned
% they are misalighned becuase of different sized windows and different
% conventions in deciding what to call the first window
% PC: 16ms clipped at start and end
% Spec/sm: 8ms clipped.
% Therefore will change PC timebins to start and end 8ms inside Spec
% range
for ll=1:length(Params.FieldsToCheck);
    NumPCtbins=length(Params.tf_bins.tPC);
    fieldname=Params.FieldsToCheck{ll};
    Params.tf_bins.tPC=linspace(Params.tf_bins.tSP(1)+0.008,Params.tf_bins.tSP(end)-0.008,NumPCtbins);
end




%% Get Stim times (i.e. to plot)
notenum=Params.notenum_stim; % which note is trig?
notefield=(['Note' num2str(notenum)]);

for ii=1:NumFields;
    fieldname=Params.FieldsToCheck{ii};
    NumSyls=length(DatStructCompiled.(fieldname));
    
    for i=1:NumSyls;
        if ~isempty(DatStructCompiled.(fieldname)(i).MostRecentTrig.(notefield).TimeSince); % since sometimes no trigger at all
            StatsStruct.(fieldname).TimeSinceLastTrig(i)=DatStructCompiled.(fieldname)(i).MostRecentTrig.(notefield).TimeSince;
        else
            StatsStruct.(fieldname).TimeSinceLastTrig(i)=nan;
        end
    end
end

%% PLOT - each field separately
tPC=Params.tf_bins.tPC;
tSP=Params.tf_bins.tSP;
fSP=Params.tf_bins.fSP;

for iii=1:NumFields;
    
    fieldname=Params.FieldsToCheck{iii};
    NumSyls=length(DatStructCompiled.(fieldname));
    figure; hold on;
    
    % 1) PC
    PC=StatsStruct.(fieldname).PC;
    
    
    h1(1)=subplot(4,1,1); hold on;
    plot(tPC*1000,PC,'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
    plot(tPC*1000,mean(PC'),'Linewidth',2)
    ylabel('Frequency (Hz)')
    title([fieldname ': Pitch contours, Specgram, W.Entropy, and Amplitude. n=' num2str(NumSyls) '.']);
    xlim([tPC(1)*1000 tPC(end)*1000])
    
    % 2) SPEC
    sp_mean=StatsStruct.(fieldname).sp_mean;
    
    
    % first, convert any sp values of 0 to non-zero(to the lowest value present);
    % solves problem of taking log of 0
    pp=find(sp_mean>0);
    mntmp = min(min(sp_mean(pp)));
    pp=find(sp_mean==0);
    sp_mean(pp) = mntmp;
    
    % second, take log
    sptemp=log(sp_mean);
    sptemp = sptemp - min(min(sptemp));
    sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
    
    h1(2)=subplot(4,1,2); hold on;
    imagesc(tSP*1000, fSP, sptemp);
    ylabel('Frequency (hz)');
    axis([tSP(1) tSP(end) fSP(1) fSP(end)]);
    
    
    
    % 3) Weiner Entropy
    h1(3)=subplot(4,1,3); hold on;
    WE=StatsStruct.(fieldname).WEntropyTimecourse;
    
    plot(tSP*1000,WE,'LineStyle','--','Color',[0.6 0.6 0.6]);
    plot(tSP*1000,mean(WE,2),'Linewidth',2,'Color','r'); % mean
    
    ylabel('(log) Weiner Entropy (-inf to 0)');
    
    
    
    % 4) Amplitude
    sm_log_mean=StatsStruct.(fieldname).sm_log_mean;
    sm_log=StatsStruct.(fieldname).sm_log;

    tSM=linspace(tSP(1),tSP(end),length(sm_log_mean));
    h1(4)=subplot(4,1,4); hold on;
    
    % individual contours
    plot(tSM*1000,sm_log,'LineStyle','--','Color',[0.6 0.6 0.6]); % individual contours
    plot(tSM*1000,sm_log_mean','Linewidth',2); % mean
    ylabel('Smoothed Amplitude (log scale)')
    xlim([tSP(1) tSP(end)])
    xlabel(['Time (ms); aligned at ' num2str(PreDur) 'ms'])
    
    
    lt_Opto_Stim_analy_PLOT_Compare2_PLOTSTIMS
    linkaxes(h1,'x');
    
end

Params.tf_bins.tSM=tSM;

%% PLOT - OVERLAY FIELDS ON SAME PLOT
tPC=Params.tf_bins.tPC;
tSP=Params.tf_bins.tSP;


if NumFields>1;
    
    PlotCols=lt_make_plot_colors(NumFields,0);
    figure; hold on;
    
    % First, PC
    for i=1:NumFields;
        fieldname=Params.FieldsToCheck{i};
        PC=StatsStruct.(fieldname).PC;
        NumSyls=size(StatsStruct.(fieldname).PC,2);
        
        % Plot all PCs
        subplot(3,1,1); hold on;
        
        plot(tPC*1000,PC,'LineStyle','-','Color',PlotCols{i},'LineWidth',0.1) % plot all pitch contours in light shade
        ylabel('Frequency (Hz)')
        title('Individual Pitch Contours');
        xlabel(['Time (ms); aligned at ' num2str(PreDur) ' ms']);
        
        xlim([tPC(1)*1000 tPC(end)*1000])
        
        
        
        % Plot mean PC (+ SEM)
        subplot(3,1,2); hold on;
        PC_sd=std(PC');
        PC_sem=PC_sd/sqrt(NumSyls-1);
        
        shadedErrorBar(tPC*1000,mean(PC'),PC_sem,{'Linewidth',2,'Color',PlotCols{i}},1);
        
        %         % annotate SEM
        %         hfig2(i)=plot(tPC*1000,mean(PC')+PC_sem,'Linewidth',0.5,'Color',PlotCols{i});
        %         plot(tPC*1000,mean(PC')-PC_sem,'Linewidth',0.5,'Color',PlotCols{i});
        
        
        xlim([tPC(1)*1000 tPC(end)*1000])
        title('Mean +/- SEM  Pitch Contours');
        ylabel('Frequency (Hz)')
        xlim([tPC(1)*1000 tPC(end)*1000])
        
        % PlOT STIMS
        lt_Opto_Stim_analy_PLOT_Compare2_PLOTSTIMS
        
        % PLOT STD
        subplot(3,1,3); hold on;
        hfig2(i)=plot(tPC*1000,PC_sd,'Linewidth',2,'Color',PlotCols{i});
        title('Pitch contour STD');
        ylabel('Pitch STD (hz)');
        xlabel(['Time (ms); aligned at ' num2str(PreDur) ' ms']);
        xlim([tPC(1)*1000 tPC(end)*1000])
        
        
        % PUT INTO STATS STRUCT
        StatsStruct.(fieldname).PC_sd=PC_sd;
        StatsStruct.(fieldname).PC_sem=PC_sem;
        
        
    end
    legend(hfig2,Params.FieldsToCheck);
    
    
    
    % Second, amplitude
    figure; hold on;
    for i=1:NumFields;
        fieldname=Params.FieldsToCheck{i};
        NumSyls=size(StatsStruct.(fieldname).PC,2);
        
        sm_log_mean=StatsStruct.(fieldname).sm_log_mean;
        sm_log=StatsStruct.(fieldname).sm_log;
        
        
        tSM=linspace(tSP(1),tSP(end),length(sm_log_mean));
        
        
        % Plot all contours
        subplot(3,1,1); hold on;
        plot(tSM*1000,sm_log,'LineStyle','-','Color',PlotCols{i},'LineWidth',0.1); % individual contours
        
        ylabel('Amplitude (log scale)')
        xlim([tSM(1) tSM(end)*1000])
        title('Individual Amplitude contours');
        
        
        % Plot Mean contours
        subplot(3,1,2); hold on;
        
        % Get STD of sm log
        sm_log_SD=std(sm_log');
        sm_log_SEM=sm_log_SD/sqrt(NumSyls-1);
        
        %     h3(i)=plot(tSM*1000,log(sm_mean'),'Linewidth',3,'Color',PlotCols{i}); % mean
        %     plot(tSM*1000,log(sm_mean')+sm_log_SD,'Linewidth',0.5,'Color',PlotCols{i}); % mean
        %     plot(tSM*1000,log(sm_mean')-sm_log_SD,'Linewidth',0.5,'Color',PlotCols{i}); % mean
        
        shadedErrorBar(tSM*1000,sm_log_mean',sm_log_SEM,{'Linewidth',2,'Color',PlotCols{i}},1); % mean
        
        ylabel('Amplitude (log scale)')
        xlim([tSM(1) tSM(end)*1000])
        title('Mean +/- SEM amplitude contours');
        
        % PLOT STIMS
        lt_Opto_Stim_analy_PLOT_Compare2_PLOTSTIMS
        
        
        % PLOT STD
        subplot(3,1,3); hold on;
        hfig3(i)=plot(tSM*1000,sm_log_SD,'Linewidth',2,'Color',PlotCols{i});
        xlabel(['Time (ms); aligned at ' num2str(PreDur) 'ms.'])
        ylabel('STD of Amplitude');
        title('STD of amplitude');
        xlim([tSM(1) tSM(end)*1000])
        
        
        % PUT INTO STATS STRUCT
        StatsStruct.(fieldname).sm_log_SD=sm_log_SD;
        StatsStruct.(fieldname).sm_log_SEM=sm_log_SEM;
        
        
    end
    
    legend(hfig3,Params.FieldsToCheck);
    
    
    % THIRD - Weiner Entropy
    figure; hold on;
    for i=1:NumFields;
        fieldname=Params.FieldsToCheck{i};
        NumSyls=size(StatsStruct.(fieldname).PC,2);
        
        WE=StatsStruct.(fieldname).WEntropyTimecourse;
        WE_mean=mean(WE,2);
        
        % Get STD
        WE_std=std(WE');
        WE_sem=WE_std./sqrt(NumSyls-1);
        
        
        % Plot all contours
        subplot(3,1,1); hold on;
        plot(tSP*1000,WE,'LineStyle','-','Color',PlotCols{i},'LineWidth',0.1); % individual contours
        
        ylabel('W. Entropy (log: -inf:0)')
        xlim([tSP(1) tSP(end)*1000])
        title('Individual W. Entropy contours');
        
        
        % Plot Mean contours
        subplot(3,1,2); hold on;
        
        shadedErrorBar(tSP*1000,WE_mean,WE_sem,{'Linewidth',2,'Color',PlotCols{i}},1); % mean
        
        
        ylabel('W. Entropy (log, -inf:0)')
        xlim([tSP(1) tSP(end)*1000])
        title('Mean +/- SEM W. Entropy contours');
        
        
        
        % PLOT STIMS
        lt_Opto_Stim_analy_PLOT_Compare2_PLOTSTIMS
        
        
        % PLOT STD
        subplot(3,1,3); hold on;
        hfig4(i)=plot(tSP*1000,WE_std,'Linewidth',2,'Color',PlotCols{i});
        xlabel(['Time (ms); aligned at ' num2str(PreDur) 'ms.'])
        xlim([tSP(1) tSP(end)*1000])
        ylabel('W. Entropy STD');
        title('STD');
        
        
        % PUT STD AND SEM INTO DATA STRUCTURE
        StatsStruct.(fieldname).WE_std=WE_std;
        StatsStruct.(fieldname).WE_sem=WE_sem;
        
    end
    legend(hfig4,Params.FieldsToCheck);
end

%% REMOVE spectrograms, takes up too much memory

for k=1:NumFields;
    fieldname=Params.FieldsToCheck{k};

    StatsStruct.(fieldname)=rmfield(StatsStruct.(fieldname),'sp');
end

%% REMOVE sm_log, save memory

for k=1:NumFields;
    fieldname=Params.FieldsToCheck{k};
    StatsStruct.(fieldname)=rmfield(StatsStruct.(fieldname),'sm_log');
end


%% SAVE;
disp('Saving...');
% go to savefolder
cd(Params.savefolder);
tstamp=lt_get_timestamp(0);

% Save StatsStruct and overwrite any old version. since these are all made
% from the same raw data.
save('StatsStruct.mat','StatsStruct');
save('Params.mat','Params');

% write down what was done in txt file name
DoneNote='DONE_Compare2';
for i=1:length(Params.FieldsToCheck);
    DoneNote=[DoneNote '_' Params.FieldsToCheck{i}];
end
DoneNote=[DoneNote '_' tstamp '.txt'];

fid1=fopen(DoneNote,'w');
fclose(fid1);

% savefigs
try
    cd('FIGURES/Compare2');
catch err
    mkdir('FIGURES/Compare2');
    cd('FIGURES/Compare2');
end

lt_save_all_figs
cd ('../../')

disp('Done!')



end


