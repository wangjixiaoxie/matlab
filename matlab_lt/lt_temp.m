clear all; close all;


% PARAMETERS
syllables={'b','c'};
NumSyls=length(syllables);
batch='batch.labeled.all';

frequency_range={[2800 4000],[2000 4000]}; % for pitch contour

[birdname bluejaynum date phrase]=lt_get_birdname_date_from_dir(1);

timestamp=lt_get_timestamp;
curr_dir=pwd;

%% FIRST, collect syllable data - sound data, pitch contour, spectrogram, timestamps, sequence

% for all single syllables

calcFFT=0; % don't get FFT at arbitrary timeshift
pc_window=[0.06,0.1]; % size of syl data window (in sec);
fs=32000;

for i=1:NumSyls;
    fvalsstr_forpc_all.(syllables{i}) =...
        findwnote2tw_v4_LT(batch, syllables{i},'',0,...
        frequency_range{i},pc_window(i)*fs,1,1,'obs0',calcFFT); % used to be 8000 but now 3000 as it is much faster in pc step (i.e. how much data to look at)
    
    NumRends(i)=length(fvalsstr_forpc_all.(syllables{i}));
    
    % PUT INTO OUTPUT STRUCTURE
    for ii=1:NumRends(i); % all renditions
        %     compiled_data.data.(syllables{i}){:,1}; % FF
        %     compiled_data.data.(syllables{i}){:,2}; % pc
        compiled_data.data.(syllables{i}){ii,3}=fvalsstr_forpc_all.(syllables{i})(ii).datt;
        %     compiled_data.data.(syllables{i}){:,4} % spectrogram
        compiled_data.data.(syllables{i}){ii,5}=fvalsstr_forpc_all.(syllables{i})(ii).fn;
        compiled_data.data.(syllables{i}){ii,6}=fvalsstr_forpc_all.(syllables{i})(ii).datenum;
        compiled_data.data.(syllables{i}){ii,7}=fvalsstr_forpc_all.(syllables{i})(ii).lbl;
        compiled_data.data.(syllables{i}){ii,8}=fvalsstr_forpc_all.(syllables{i})(ii).NotePos;
        compiled_data.data.(syllables{i}){ii,9}=fvalsstr_forpc_all.(syllables{i})(ii).TRIG;
    end
end



%% Second, COLLECT PITCH CONTOUR DATA

pc_harms=1;
pc_time_window={[25 120],[45 220]};
plot_result=1;

for i=1:NumSyls;
    [PC_raw.(syllables{i}), pc_F(i), pc_T(i), spec]=jc_pitchcontourFV_LT(fvalsstr_forpc_all.(syllables{i}),...
        1024,1020,1, frequency_range{i}(1),...
        frequency_range{i}(2),pc_harms,'obs0');
    
    for ii=1:NumRends(i); % all renditions
        compiled_data.data.(syllables{i}){ii,2}=PC_raw.(syllables{i})(:,ii)'; % pc
        compiled_data.data.(syllables{i}){ii,4}=spec{ii}; % spectrogram
    end
end

%% PLOT PITCH CONTOUR
if plot_result==1;
    for i = 1:NumSyls;
        %plots the mean and individual renditions
        figure; hold on
        plot(PC_raw.(syllables{i}),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(mean(PC_raw.(syllables{i})'),'Linewidth',2)
        
        xlabel(['Time, bins: ' num2str(pc_T(2)-pc_T(1))])
        ylabel('Frequency (Hz)')
        title(['Mean and individual pitch contours for ' syllables{i}]);
        
        % plot with real time axis
        figure; hold on
        plot(pc_T(i)*1000,PC_raw.(syllables{i}),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(pc_T(i)*1000,mean(PC_raw.(syllables{i})'),'Linewidth',2)
        xlabel('Time (ms)')
        ylabel('Frequency (Hz)')
        title(['Mean and individual pitch contours for ' syllables{i}]);
    end
end


%% REMOVING OUTLIERS (i.e. pc that fluctuates wildly past threshold)
% have to know what time bins you want to use to calculate pitch. for each
% time bin, will perform this outlier test.

for i=1:length(syllables);
    PC_mat=cell2mat(compiled_data.data.(syllables{i})(:,2)); % convert to mat
    columns_to_check=pc_time_window{i};
    num_columns=columns_to_check(2)-columns_to_check(1)+1; % how many times to do analysis
    Outliers{i}=[];
    
    for ii=1:num_columns; % for each column, perform check.
        [~, B, C]=lt_db_tukey_outlier(PC_mat,columns_to_check(1)+ii-1,2.5);
        Outliers{i}=[Outliers{i} B C]; % accumulate index of outliers
    end
    OutlInds{i}=unique(Outliers{i}); % gives back specific index only once
end

% Remove outliers from data:
compiled_data.data_NoOutlier=compiled_data.data;
for i=1:length(syllables);
    compiled_data.data_NoOutlier.(syllables{i})(OutlInds{i},:)=[]; % remove outlier rows
    disp(['Removed ' num2str(length(OutlInds{i})) ' outlier (based on PC) for syllable ' syllables{i} ]);
end


% PLOT TO CONFIRM THAT OUTLIERS REMOVED:
for i = 1:length(syllables);
    if ~isempty(OutlInds{i});
        
        %plots the mean and individual renditions
        figure; hold on
        PCmat=cell2mat(compiled_data.data_NoOutlier.(syllables{i})(:,2));
        plot(PCmat','LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(mean(PC_raw.(syllables{i})'),'Linewidth',2)
        xlabel(['Time, bins: ' num2str(pc_T(2)-pc_T(1))])
        ylabel('Frequency (Hz)')
        title(['OUTLIER GONE: Mean and individual pitch contours for ' syllables{i}]);
        
        % PLOTS THE OUTLIERS ONLY
        %plots the mean and individual renditions
        figure; hold on
        PCmat=cell2mat(compiled_data.data.(syllables{i})(OutlInds{i},2));
        plot(PCmat','LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        xlabel(['Time, bins: ' num2str(pc_T(2)-pc_T(1))])
        ylabel('Frequency (Hz)')
        title(['OUTLIERS THAT WERE REMOVED: ' syllables{i}]);
    end
end

%% EXTRACT FF DATA FROM PITCH CONTOUR


% ALL DATA (including outliers)
for i=1:length(syllables);
    PCmat=cell2mat(compiled_data.data.(syllables{i})(:,2)); % cell to mat
    PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
    
    N=NumRends(i); % num rends
    compiled_data.data.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
end

% OUTLIERS REMOVED
for i=1:length(syllables);
    PCmat=cell2mat(compiled_data.data_NoOutlier.(syllables{i})(:,2)); % cell to mat
    PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
    
    N=size(compiled_data.data_NoOutlier.(syllables{i}),1); % num rends
    compiled_data.data_NoOutlier.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
end


%% FILTER syllables based on sequence context (e.g. only b from ab);

% Initiate with list of sequences you want to extract
SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should combine
SylTargList={'c','b'}; % these must already have raw data compiled above
SeqPostList={'b',''};

for j=1:length(SylTargList);
    SeqPre=SeqPreList{j};
    SylTarg=SylTargList{j};
    SeqPost=SeqPostList{j};
    
    % Automatic params
    TargStr=[SeqPre SylTarg SeqPost]; % entire string
    TartPos=length(SeqPre)+1; % target syl position in string
    
    
    % OUTLIERS PRESENT
    NumRends=size(compiled_data.data.(SylTarg),1);
    IndsToKeep=[];
    for ii=1:NumRends; % all rends of the targ syl
        sylpos=compiled_data.data.(SylTarg){ii,8}; % ind in song of candidate
        IndsToCheck=(sylpos-length(SeqPre)):(sylpos+length(SeqPost)); % what within song inds to look at
        SeqToCheck=compiled_data.data.(SylTarg){ii,7}(IndsToCheck); % is this the desired string?
        
        if strcmp(TargStr,SeqToCheck)==1; % if is desired string, then keep this datapt
            IndsToKeep=[IndsToKeep ii];
        end
    end
    
    compiled_data.data_SeqFiltered.([SeqPre '_' SylTarg '_' SeqPost])=...
        compiled_data.data.(SylTarg)(IndsToKeep,:);
    
    
    
    % OUTLIERS REMOVED
    NumRends=size(compiled_data.data_NoOutlier.(SylTarg),1);
    IndsToKeep=[];
    for ii=1:NumRends;
        sylpos=compiled_data.data_NoOutlier.(SylTarg){ii,8}; % ind in song of candidate
        IndsToCheck=(sylpos-length(SeqPre)):(sylpos+length(SeqPost)); % what within song inds to look at
        SeqToCheck=compiled_data.data_NoOutlier.(SylTarg){ii,7}(IndsToCheck); % is this the desired string?
        
        if strcmp(TargStr,SeqToCheck)==1; % if is desired string, then keep this datapt
            IndsToKeep=[IndsToKeep ii];
        end
    end
    
    compiled_data.data_SeqFiltered_NoOutliers.([SeqPre '_' SylTarg '_' SeqPost])=...
        compiled_data.data_NoOutlier.(SylTarg)(IndsToKeep,:);
end



%% GET SUMMARY STATS ACROSS ALL RENDS

% RAW (outliers present)
for i=1:NumSyls;
    X=cell2mat(compiled_data.data.(syllables{i})(:,1)); % FF vals to matrix
    compiled_data.summary_stats.(syllables{i}).meanFF=mean(X);
    compiled_data.summary_stats.(syllables{i}).medianFF=median(X);
    compiled_data.summary_stats.(syllables{i}).sdFF=std(X);
    compiled_data.summary_stats.(syllables{i}).n=length(X);
    compiled_data.summary_stats.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_data.data.(syllables{i})(:,2)); % rows: rend, cols: timebin
    compiled_data.summary_stats.(syllables{i}).pitchcountour_mean=mean(X,1);
    
    % get mean spectrogram
    X=[];
    for ii=1:size(compiled_data.data.(syllables{i})(:,4),1); % num rends
        X(:,:,ii)=compiled_data.data.(syllables{i}){ii,4}; % extract individual specs
    end
    compiled_data.summary_stats.(syllables{i}).spec_mean=mean(X,3);
end

% OUTLIER REMOVED
for i=1:NumSyls;
    X=cell2mat(compiled_data.data_NoOutlier.(syllables{i})(:,1)); % FF vals to matrix
    compiled_data.summary_stats_NoOutlier.(syllables{i}).meanFF=mean(X);
    compiled_data.summary_stats_NoOutlier.(syllables{i}).medianFF=median(X);
    compiled_data.summary_stats_NoOutlier.(syllables{i}).sdFF=std(X);
    compiled_data.summary_stats_NoOutlier.(syllables{i}).n=length(X);
    compiled_data.summary_stats_NoOutlier.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_data.data_NoOutlier.(syllables{i})(:,2)); % rows: rend, cols: timebin
    compiled_data.summary_stats_NoOutlier.(syllables{i}).pitchcountour_mean=mean(X,1);
    
    % get mean spectrogram
    X=[];
    for ii=1:size(compiled_data.data_NoOutlier.(syllables{i})(:,4),1); % num rends
        X(:,:,ii)=compiled_data.data_NoOutlier.(syllables{i}){ii,4}; % extract individual specs
    end
    compiled_data.summary_stats_NoOutlier.(syllables{i}).spec_mean=mean(X,3);

end




% PERFORM SAME FOR SEQ-DEP DATA
syl_fields=fields(compiled_data.data_SeqFiltered); % what seqs are there?
for i=1:length(syl_fields);

    % OUTLIERS PRESENT
    X=cell2mat(compiled_data.data_SeqFiltered.(syl_fields{i})(:,1));
    compiled_data.summary_stats.(syl_fields{i}).meanFF=mean(X);
    compiled_data.summary_stats.(syl_fields{i}).medianFF=median(X);
    compiled_data.summary_stats.(syl_fields{i}).sdFF=std(X);
    compiled_data.summary_stats.(syl_fields{i}).n=length(X);
    compiled_data.summary_stats.(syl_fields{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_data.data_SeqFiltered.(syl_fields{i})(:,2)); % rows: rend, cols: timebin
    compiled_data.summary_stats.(syl_fields{i}).pitchcountour_mean=mean(X,1);

    % get mean spectrogram
    X=[];
    for ii=1:size(compiled_data.data_SeqFiltered.(syl_fields{i})(:,4),1); % num rends
        X(:,:,ii)=compiled_data.data_SeqFiltered.(syl_fields{i}){ii,4}; % extract individual specs
    end
    compiled_data.summary_stats.(syl_fields{i}).spec_mean=mean(X,3);

    
    
    % NO OUTLIERS
    X=cell2mat(compiled_data.data_SeqFiltered_NoOutliers.(syl_fields{i})(:,1));
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).meanFF=mean(X);
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).medianFF=median(X);
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).sdFF=std(X);
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).n=length(X);
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_data.data_SeqFiltered_NoOutliers.(syl_fields{i})(:,2)); % rows: rend, cols: timebin
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).pitchcountour_mean=mean(X,1);
    
    % get mean spectrogram
    X=[];
    for ii=1:size(compiled_data.data_SeqFiltered_NoOutliers.(syl_fields{i})(:,4),1); % num rends
        X(:,:,ii)=compiled_data.data_SeqFiltered_NoOutliers.(syl_fields{i}){ii,4}; % extract individual specs
    end
    compiled_data.summary_stats_NoOutlier.(syl_fields{i}).spec_mean=mean(X,3);
end



%% SAVE DATA

savedir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/


%% troubleshooting

% TO FIND BINS IN PC THAT ARE PAST A CERTAIN THRESHOLD.
% for i=1:size(PC_raw.b,2);
%     if sum(find(PC_raw.b(1:150,i)>3600))>1;
%         disp(i);
%         disp(find(PC_raw.b(1:150,i)>3600))
%     end
% end
%
