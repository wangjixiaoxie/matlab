function compiled_seqdep_pitch=lt_compile_seq_dep_pitch_data(batch, syllables, frequency_range, pc_window, pc_time_window, plotON)
%% temp: removing things to make file smaller - i.e. dont save spectrograms (just average).
%% LT - 11/13/14 -
% Compiles pitch data - pitch contour, summary pitch data, and sequence
% data. Also
% removes outliers based on the pitch contour (i.e. any timebin with
% extreme values removes the entire syllable data pt).

% Once this is run, you can use lt_compile_seq_dep_pitch_data_SEQFILTER to
% filter the data in a sequence dependent
% manner (e.g. want not just b, but b that follows "a" specifically.
% Run in day folder.  Makes structure and saves in a folder one dir up.

% examples inputs:
% batch='batch.labeled.all';
% syllables={'b','c'};
% frequency_range={[3000 3900],[2000 3000]}; % for findwnote
% pc_window=[0.06,0.1]; % size of syl data window (in sec); (how much data to get for each rend (relative to onset), in sec)
% pc_time_window={[25 120],[45 220]}; for pitch contour (time bins to avg).
% SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should combine
% SylTargList={'c','b'}; % these must already have raw data compiled above
% SeqPostList={'b',''};
% plotON=1;


% PARAMETERS
% findwnote
calcFFT=0; % don't get FFT at arbitrary timeshift
fs=32000;

% pitch contour
pc_harms=1; % harmonics to take weighted avg over. 1 or 2 is good.
plotON=1;

% AUTO params
NumSyls=length(syllables);
[birdname bluejaynum date phrase]=lt_get_birdname_date_from_dir(1);
timestamp=lt_get_timestamp;
curr_dir=pwd;


%% FIRST, collect syllable data - sound data, pitch contour, spectrogram, timestamps, sequence
SylsToDrop=[]; % will fill up with syls to skip (no data);
NumSyls2=NumSyls;
for i=1:NumSyls;
    try % i.e. maybe this day does not have this syllable
        fvalsstr_forpc_all.(syllables{i}) =...
            findwnote2tw_v4_LT(batch, syllables{i},'',-0.005,...
            frequency_range{i},pc_window(i)*fs,1,1,'obs0',calcFFT);
        
        NN=length(fvalsstr_forpc_all.(syllables{i}));
        
        % PUT INTO OUTPUT STRUCTURE
        for ii=1:NN; % all renditions
            compiled_seqdep_pitch.data.(syllables{i}){ii,3}=fvalsstr_forpc_all.(syllables{i})(ii).datt;
            compiled_seqdep_pitch.data.(syllables{i}){ii,5}=fvalsstr_forpc_all.(syllables{i})(ii).fn;
            compiled_seqdep_pitch.data.(syllables{i}){ii,6}=fvalsstr_forpc_all.(syllables{i})(ii).datenum;
            compiled_seqdep_pitch.data.(syllables{i}){ii,7}=fvalsstr_forpc_all.(syllables{i})(ii).lbl;
            compiled_seqdep_pitch.data.(syllables{i}){ii,8}=fvalsstr_forpc_all.(syllables{i})(ii).NotePos;
            compiled_seqdep_pitch.data.(syllables{i}){ii,9}=fvalsstr_forpc_all.(syllables{i})(ii).TRIG;
        end
        
    catch err % if syl does not exist, then remove it from analysis.
        NumSyls2=NumSyls2-1;
        SylsToDrop=[SylsToDrop i];
    end
    
end

syllables(SylsToDrop)=[]; % removing NoData syls.
frequency_range(SylsToDrop)=[];
pc_window(SylsToDrop)=[];
pc_time_window(SylsToDrop)=[];

NumSyls=NumSyls2;

% how many renditions for each syl.
for i=1:NumSyls;
    NumRends(i)=length(fvalsstr_forpc_all.(syllables{i}));
end


%% Second, COLLECT PITCH CONTOUR DATA

for i=1:NumSyls;
    [PC_raw.(syllables{i}), pc_F{i}, pc_T{i}, ~]=jc_pitchcontourFV_LT(fvalsstr_forpc_all.(syllables{i}),...
        1024,1020,1, frequency_range{i}(1),...
        frequency_range{i}(2),pc_harms,'obs0');
    
    for ii=1:NumRends(i); % all renditions
        compiled_seqdep_pitch.data.(syllables{i}){ii,2}=PC_raw.(syllables{i})(:,ii)'; % pc
        % removed spectrogram becuase takes too much space - for 300 syls,
        % about 30MB
        %         compiled_seqdep_pitch.data.(syllables{i}){ii,4}=spec{ii}; % spectrogram
        
    end
end


%% PLOT PITCH CONTOUR
if plotON==1;
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
    end
end


%% REMOVING OUTLIERS (i.e. pc that fluctuates wildly past threshold) - saves to separate field.
% have to know what time bins you want to use to calculate pitch. for each
% time bin, will perform this outlier test.
% WILL SAVE INDICES OF 1) OUTLIERS AND 2) NON-OUTLIERS, BUT NOT ACTUALLY
% FILTER THE DATA.

for i=1:length(syllables);
    PC_mat=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,2)); % convert to mat
    columns_to_check=pc_time_window{i};
    num_columns=columns_to_check(2)-columns_to_check(1)+1; % how many times to do analysis
    Outliers{i}=[];
    
    for ii=1:num_columns; % for each column, perform check.
        [~, B, C]=lt_db_tukey_outlier_90tile(PC_mat,columns_to_check(1)+ii-1,1.5);
        Outliers{i}=[Outliers{i} B' C']; % accumulate index of outliers
    end
    OutlInds{i}=unique(Outliers{i}); % gives back specific index only once
end


% Keep a copy of raw data (with outliers), and outlier indices
compiled_seqdep_pitch.data_WithOutlier=compiled_seqdep_pitch.data;

for i=1:length(syllables);
    compiled_seqdep_pitch.PARAMETERS.OutlierInds.(syllables{i})=OutlInds{i}; % keep indices
end


% REMOVE outliers from data:
for i=1:length(syllables);
    compiled_seqdep_pitch.data.(syllables{i})(OutlInds{i},:)=[]; % remove outlier rows
    disp(['Removed ' num2str(length(OutlInds{i})) ' outlier (based on PC) for syllable ' syllables{i} ]);
end


% PLOT TO CONFIRM THAT OUTLIERS REMOVED:
for i = 1:length(syllables);
    if ~isempty(OutlInds{i});
        
        %plots the mean and individual renditions
        figure; hold on, subplot(2,1,1), hold on;
        PCmat=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,2));
        plot(PCmat','LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(mean(PCmat),'Linewidth',2)
        xlabel(['Time, bins: ' num2str(pc_T{i}(2)-pc_T{i}(1))])
        ylabel('Frequency (Hz)')
        title(['OUTLIER GONE: Mean and individual pitch contours for ' syllables{i}]);
        
        % PLOTS THE OUTLIERS ONLY
        %plots the mean and individual renditions
        subplot(2,1,2); hold on;
        PCmat=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(OutlInds{i},2));
        plot(PCmat','LineStyle','--')
        xlabel(['Time, bins: ' num2str(pc_T{i}(2)-pc_T{i}(1))])
        ylabel('Frequency (Hz)')
        title(['OUTLIERS THAT WERE REMOVED: ' syllables{i}]);
    end
end

%% EXTRACT FF DATA FROM PITCH CONTOUR

% OUTLIERS REMOVED
for i=1:length(syllables);
    PCmat=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,2)); % cell to mat
    PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
    
    N=size(compiled_seqdep_pitch.data.(syllables{i}),1); % num rends
    compiled_seqdep_pitch.data.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
end

% WITH OUTLIERS PRESENT
for i=1:length(syllables);
    PCmat=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,2)); % cell to mat
    PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
    
    N=size(compiled_seqdep_pitch.data_WithOutlier.(syllables{i}),1); % num rends
    compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
end


%% GET SUMMARY STATS ACROSS ALL RENDS

% DATA WITHOUT OUTLIERS
for i=1:NumSyls;
    X=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,1)); % FF vals to matrix
    compiled_seqdep_pitch.summary_stats.(syllables{i}).meanFF=mean(X);
    compiled_seqdep_pitch.summary_stats.(syllables{i}).medianFF=median(X);
    compiled_seqdep_pitch.summary_stats.(syllables{i}).sdFF=std(X);
    compiled_seqdep_pitch.summary_stats.(syllables{i}).n=length(X);
    compiled_seqdep_pitch.summary_stats.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,2)); % rows: rend, cols: timebin
    compiled_seqdep_pitch.summary_stats.(syllables{i}).pitchcountour_mean=mean(X,1);
    
    %     % get mean spectrogram
    %     X=[];
    %     for ii=1:size(compiled_seqdep_pitch.data.(syllables{i})(:,4),1); % num rends
    %         if
    %         X(:,:,ii)=spec.(syllables{i}){ii}; % extract individual specs
    %     end
    %     compiled_seqdep_pitch.summary_stats.(syllables{i}).spec_mean=mean(X,3);
end

% OUTLIER REMOVED
for i=1:NumSyls;
    X=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,1)); % FF vals to matrix
    compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).meanFF=mean(X);
    compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).medianFF=median(X);
    compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).sdFF=std(X);
    compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).n=length(X);
    compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,2)); % rows: rend, cols: timebin
    compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).pitchcountour_mean=mean(X,1);
    
    % get mean spectrogram
    %     X=[];
    %     for ii=1:size(compiled_seqdep_pitch.data_NoOutlier.(syllables{i})(:,4),1); % num rends
    %         X(:,:,ii)=compiled_seqdep_pitch.data_NoOutlier.(syllables{i}){ii,4}; % extract individual specs
    %     end
    %     compiled_seqdep_pitch.summary_stats_NoOutlier.(syllables{i}).spec_mean=mean(X,3);
    
end




%% SAVE DATA

% Compile parameters
% syl specific
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.syllables=syllables;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.batch=batch;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.frequency_range=frequency_range;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.pc_window=pc_window;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.pc_F=pc_F;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.pc_T=pc_T;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.pc_harms=pc_harms;
compiled_seqdep_pitch.PARAMETERS.SINGLESYLS.pc_time_window=pc_time_window;

% global
compiled_seqdep_pitch.PARAMETERS.birdname=birdname;
compiled_seqdep_pitch.PARAMETERS.bluejaynum=bluejaynum;
compiled_seqdep_pitch.PARAMETERS.date=date{2};
compiled_seqdep_pitch.PARAMETERS.phrase=phrase;
compiled_seqdep_pitch.PARAMETERS.timestamp=timestamp;

% already assigned above
compiled_seqdep_pitch.PARAMETERS.OutlierInds;

% Assign legend (what do data columns mean?)
compiled_seqdep_pitch.PARAMETERS.Legend_data{1,:}={'1_FF','2_PitchContour','3_SoundDat','4_Spec',...
    '5_FileName','6_DateNum','7_LabelStr','8_NotePos','9_Trig'};


% SAVE
savedir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/compile_seq_dep_pitch_data_' phrase];
savename=['DataStruct' '_' date{2}];

try
    cd(savedir);
catch err % dir does not exist
    mkdir(savedir);
    cd(savedir);
end

save(savename,'compiled_seqdep_pitch','-v7.3');

cd(curr_dir)

disp('DONE!');

%% troubleshooting

% TO FIND BINS IN PC THAT ARE PAST A CERTAIN THRESHOLD.
% for i=1:size(PC_raw.b,2);
%     if sum(find(PC_raw.b(1:150,i)>3600))>1;
%         disp(i);
%         disp(find(PC_raw.b(1:150,i)>3600))
%     end
% end
%
