function Params=lt_seq_dep_pitch_SeqFilter(Params)

%% LT 1/31/15 - compiles all days into one structure, and removes outliers.

% RUN IN FOLDER CONTAINING THE saved day data structures

% This program opens all files and puts in order into one structure (based
% on date).

% input1=1; % will grab all files that exist and start structure at
% earliest day

% or:

% input1=0; % will have to specifify first and last day
% FirstDay='11Oct2014';
% LastDay='02Nov2014';

%%
% COLLECT ALL DATA into batch
batch=lt_write_all_folder_contents_to_batch(0);


% GO LINE BY LINE AND GET DATE AND FILE NAME
fid=fopen(batch);
c=1;
nextline=fgetl(fid);

while nextline~=-1; % i.e while files still exist
    filename{c}=nextline;
    filedate{c}=nextline(14:22);
    
    nextline=fgetl(fid);
    c=c+1;
end


% WHAT IS INDEX OF FIRST AND LAST DAY?
if input1==0; % then will use specified dates
    
    FirstDateNum=datenum(FirstDay,'ddmmmyyyy');
    LastDateNum=datenum(LastDay,'ddmmmyyyy');
elseif input1==1; % will figure out all dates
    
    FirstDateNum=min(datenum(filedate,'ddmmmyyyy'));
    LastDateNum=max(datenum(filedate,'ddmmmyyyy'));
    
    FirstDay=datestr(FirstDateNum,'ddmmmyyyy');
    LastDay=datestr(LastDateNum,'ddmmmyyyy');
end

% OPEN EACH FILE IN AND PLACE IN CORRECT POSITION (by date) IN LARGER
% STRUCTURE
AllDays_compiled_seqdep_pitch=cell(LastDateNum-FirstDateNum+1,1); % one cell for every desired day.
for i=1:length(filename);
    X=load(filename{i});
    
    ii=datenum(filedate{i},'ddmmmyyyy')-FirstDateNum+1; % date position
    
    if ii>0 && ii<=LastDateNum-FirstDateNum+1; % only compile data if it is within range of desired dates (only matters if I specified input dates)
    AllDays_compiled_seqdep_pitch{ii}=X.compiled_seqdep_pitch; % slot data into correct position
    end
    
end


%% SAVE

SaveName=['AllDays_Compiled_' FirstDay '_to_' LastDay];

try
    cd('AllDays_Compiled');
catch err
    mkdir('AllDays_Compiled');
    cd('AllDays_Compiled');
end

save(SaveName,'AllDays_compiled_seqdep_pitch','-v7.3');

cd ..

disp('DONE!');





%% outliers



% %% REMOVING OUTLIERS (i.e. pc that fluctuates wildly past threshold) - saves to separate field.
% % have to know what time bins you want to use to calculate pitch. for each
% % time bin, will perform this outlier test.
% % WILL SAVE INDICES OF 1) OUTLIERS AND 2) NON-OUTLIERS, BUT NOT ACTUALLY
% % FILTER THE DATA.
% 
% for i=1:length(syllables);
%     PC_mat=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,2)); % convert to mat
%     columns_to_check=pc_time_window{i};
%     num_columns=columns_to_check(2)-columns_to_check(1)+1; % how many times to do analysis
%     Outliers{i}=[];
%     
%     for ii=1:num_columns; % for each column, perform check.
%         [~, B, C]=lt_db_tukey_outlier_90tile(PC_mat,columns_to_check(1)+ii-1,1.5);
%         Outliers{i}=[Outliers{i} B' C']; % accumulate index of outliers
%     end
%     OutlInds{i}=unique(Outliers{i}); % gives back specific index only once
% end
% 
% 
% % Keep a copy of raw data (with outliers), and outlier indices
% compiled_seqdep_pitch.data_WithOutlier=compiled_seqdep_pitch.data;
% 
% for i=1:length(syllables);
%     compiled_seqdep_pitch.PARAMETERS.OutlierInds.(syllables{i})=OutlInds{i}; % keep indices
% end
% 
% 
% % REMOVE outliers from data:
% for i=1:length(syllables);
%     compiled_seqdep_pitch.data.(syllables{i})(OutlInds{i},:)=[]; % remove outlier rows
%     disp(['Removed ' num2str(length(OutlInds{i})) ' outlier (based on PC) for syllable ' syllables{i} ]);
% end
% 
% 
% % PLOT TO CONFIRM THAT OUTLIERS REMOVED:
% for i = 1:length(syllables);
%     if ~isempty(OutlInds{i});
%         
%         %plots the mean and individual renditions
%         figure; hold on, subplot(2,1,1), hold on;
%         PCmat=cell2mat(compiled_seqdep_pitch.data.(syllables{i})(:,2));
%         plot(PCmat','LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
%         plot(mean(PCmat),'Linewidth',2)
%         xlabel(['Time, bins: ' num2str(pc_T{i}(2)-pc_T{i}(1))])
%         ylabel('Frequency (Hz)')
%         title(['OUTLIER GONE: Mean and individual pitch contours for ' syllables{i}]);
%         
%         % PLOTS THE OUTLIERS ONLY
%         %plots the mean and individual renditions
%         subplot(2,1,2); hold on;
%         PCmat=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(OutlInds{i},2));
%         plot(PCmat','LineStyle','--')
%         xlabel(['Time, bins: ' num2str(pc_T{i}(2)-pc_T{i}(1))])
%         ylabel('Frequency (Hz)')
%         title(['OUTLIERS THAT WERE REMOVED: ' syllables{i}]);
%     end
% end



%% want to output fields like this;

% for i=1:length(syllables);
%     PCmat=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,2)); % cell to mat
%     PCff=mean(PCmat(:,pc_time_window{i}(1):pc_time_window{i}(2)),2); % take mean across time
%     
%     N=size(compiled_seqdep_pitch.data_WithOutlier.(syllables{i}),1); % num rends
%     compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
% end
% 
% 


% OUTLIER REMOVED
% for i=1:NumSyls;
%     X=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,1)); % FF vals to matrix
%     compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).meanFF=mean(X);
%     compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).medianFF=median(X);
%     compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).sdFF=std(X);
%     compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).n=length(X);
%     compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
%     
%     % get mean pitch contour
%     X=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syllables{i})(:,2)); % rows: rend, cols: timebin
%     compiled_seqdep_pitch.summary_stats_WithOutlier.(syllables{i}).pitchcountour_mean=mean(X,1);
%     
%     % get mean spectrogram
%     %     X=[];
%     %     for ii=1:size(compiled_seqdep_pitch.data_NoOutlier.(syllables{i})(:,4),1); % num rends
%     %         X(:,:,ii)=compiled_seqdep_pitch.data_NoOutlier.(syllables{i}){ii,4}; % extract individual specs
%     %     end
%     %     compiled_seqdep_pitch.summary_stats_NoOutlier.(syllables{i}).spec_mean=mean(X,3);
%     
% end
% 

