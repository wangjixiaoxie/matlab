function [Params, AllDays_RawDatStruct]=lt_seq_dep_pitch_SeqFilterCompile(Params, plotON)
%% LT 4/28/15 - 
% moved old version to '_OLD'. Here changing saving directory format to
% save structures all in one folder and also overwrite.



%% LT 3/22/15 - added ability to give syl and designate it as repeats and function will get max repeats for day
% e.g. Repeats = {'jB'}; means will look for B repeats, only after j.


%% LT 1/31/15 - modified from lt_compile_seq_dep_pitch_data_SEQFILTER_MULTDAYS
% RUN IN seq_dep_pitch folder



%% LT 11/18/14 - Loads data structures containing single syls (e.g. a, b, c), and filters out sequence specific stuff (e.g. a[b])
% Run in the folder containing the data structures (e.g. % Run in folder:
% e.g.
% /bluejay3/lucas/birds/pu11wh87/compile_seq_dep_pitch_data_SeqDepPitchShift,
% containing structures:
% e.g. DataStruct_02Nov2014);

%% PARAMETERS
% all_days = 1 (Run on all data in folder); =0 (Use the dates specified
% below). If 1, then doesn't matter what I enter for days argumemtns.
% FirstDay='11Oct2014'; %to filter out days I don't want
% LastDay='02Nov2014';

% WHAT SEQUENCES DO I WANT?
% SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should
% combine. in this case the two sequences are ac[c]b and ab[b].
% SylTargList={'c','b'}; % these must already have raw data compiled above
% SeqPostList={'b',''};


%% simplify params names
disp('Running...')

all_days=Params.SeqFilter.all_daysON; % If 1, then doesn't matter what I enter for days argumemtns.
FirstDay=Params.SeqFilter.FirstDay;
LastDay=Params.SeqFilter.LastDay;

SeqPreList=Params.SeqFilter.SeqPreList;
SylTargList=Params.SeqFilter.SylTargList;
SeqPostList=Params.SeqFilter.SeqPostList;


% Load any single day's params to get useful info

params_batch=lt_write_all_folder_contents_to_batch_v2('Params*');

ftmp=fopen(params_batch);

fname=fgetl(ftmp);
Params_singleday=load(fname);

fclose(ftmp)

% Combine old params and new params
Tmp = [fieldnames(Params)' fieldnames(Params_singleday.Params)'; struct2cell(Params)' struct2cell(Params_singleday.Params)'];
clear Params
Params=struct(Tmp{:});





%% LOAD SINGLE SYL DATA STRUCTURES

% COLLECT ALL DATA into batch
batch=lt_write_all_folder_contents_to_batch_v2('RawDatStruct*');


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
if all_days==0; % then will use specified dates
    
    FirstDateNum=datenum(FirstDay,'ddmmmyyyy');
    LastDateNum=datenum(LastDay,'ddmmmyyyy');
elseif all_days==1; % will figure out all dates
    
    FirstDateNum=min(datenum(filedate,'ddmmmyyyy'));
    LastDateNum=max(datenum(filedate,'ddmmmyyyy'));
    
    FirstDay=datestr(FirstDateNum,'ddmmmyyyy');
    LastDay=datestr(LastDateNum,'ddmmmyyyy');
end

Params.SeqFilter.FirstDay=FirstDay;
Params.SeqFilter.LastDay=LastDay;



%% PLACE EACH FILE IN CORRECT POSITION AND FILTER SEQUENCE
% STRUCTURE
AllDays_RawDatStruct=cell(LastDateNum-FirstDateNum+1,1); % one cell for every desired day.
for i=1:length(filename);
    ii=datenum(filedate{i},'ddmmmyyyy')-FirstDateNum+1; % date position (to check if I want this structure)
    
    if ii>0 && ii<=LastDateNum-FirstDateNum+1; % only compile data if it is within range of desired dates (only matters if I specified input dates)
        load(filename{i});
        [Params,RawDatStruct]=lt_seq_dep_pitch_SeqFilter_singleday(Params,RawDatStruct); % Run function on the data structure just loaded - will save to a specific folder: savedir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/compile_seq_dep_pitch_data_' phrase '/SEQFILTER'];
        
        % I used to save (in "IOLD" version of function) - not anymore as
        % never used. just use compiled structure.
%         % save 
%         savename=[Params.SeqFilter.savedir '/RawDatStruct_' filedate{i}];
%         save(savename,'RawDatStruct');
        
        
        % compile structure
        AllDays_RawDatStruct{ii}=RawDatStruct; % put dat struct in correct place
    end
end





%% REMOVE OUTLIERS (based on PC)

% have to know what time bins you want to use to calculate pitch. for each
% time bin, will perform this outlier test.
% WILL SAVE INDICES OF 1) OUTLIERS AND 2) NON-OUTLIERS, BUT NOT ACTUALLY
% FILTER THE DATA.

NumDays=datenum(Params.SeqFilter.LastDay)-datenum(Params.SeqFilter.FirstDay)+1;


% 1) for seq dep syls (e.g. ab, cb) must know what single syl it is (e.g.
% b) to know what time bins to look at. Here: figure out what syl it is.
for j=1:NumDays;
    
    % check if day has data
    if ~isempty(AllDays_RawDatStruct{j})
        syllables=fieldnames(AllDays_RawDatStruct{j}.data);

        for jj=1:length(syllables);
            
            % 1) for each seq dep syl, figure out which original syl it was
            if length(syllables{jj})>1; % only do if this is a seq.
                % go through the original syls and find the match
                actual_nt=regexp(syllables{jj},'[A-Z]','match'); % what note is the capital letter?
                actual_nt=lower(actual_nt);
                for jjj=1:length(Params.DayRawDat.syllables); % compare that syl to all possible single syls
                    if strcmp(Params.DayRawDat.syllables{jjj},actual_nt)==1;
                        Params.SeqFilter.OrigNoteID{j}(jj)=jjj; % for this day, ID all the syls
                    end
                end
            else
                actual_nt=syllables{jj}; % if syl is already single syl
                actual_nt=lower(actual_nt);
                for jjj=1:length(Params.DayRawDat.syllables);
                    if strcmp(Params.DayRawDat.syllables{jjj},actual_nt)==1;
                        Params.SeqFilter.OrigNoteID{j}(jj)=jjj; % for this day, ID all the syls
                    end
                end
                
            end
            
            % 2) Then figure out which time bins it should use
            
            noteID=Params.SeqFilter.OrigNoteID{j}(jj); % which note is it?
            Params.SeqFilter.pc_time_window_list{j}(:,jj)=Params.DayRawDat.pc_time_window{noteID}; % use that notes time windows
            
        end
        
    end
    
end



% START LOOKING FOR OUTLIERS
for j=1:NumDays;
    
    if ~isempty(AllDays_RawDatStruct{j}); % only continue if day has data.
        syllables=fieldnames(AllDays_RawDatStruct{j}.data); % syllables in this day
        Outliers={};
        OutlInds={};
        for i=1:length(syllables);
            
            PC_mat=cell2mat(AllDays_RawDatStruct{j}.data.(syllables{i})(:,2)); % convert to mat
%             PC_mat=cell2mat(AllDays_RawDatStruct{j}.data_WithOutlier.(syllables{i})(:,2)); % convert to mat

            columns_to_check=Params.SeqFilter.pc_time_window_list{j}(:,i);
            num_columns=columns_to_check(2)-columns_to_check(1)+1; % how many times to do analysis
            Outliers{i}=[];
            
            for ii=1:num_columns; % for each column, perform check.
                % empirically determined 2 to be a good iqr value to use.
                % (based on pu37wh20)
                % changed to 2.5 (6/15/15)
                [~, B, C]=lt_db_tukey_outlier_90tile(PC_mat,columns_to_check(1)+ii-1,2.5);
%                 [~, B, C]=lt_db_tukey_outlier(PC_mat,columns_to_check(1)+ii-1,1.5);

                Outliers{i}=[Outliers{i} B' C']; % accumulate index of outliers
            end
            OutlInds{i}=unique(Outliers{i}); % gives back specific index only once
        end
        
        
        
        % Remove those outliers
        AllDays_RawDatStruct{j}.data_WithOutlier=AllDays_RawDatStruct{j}.data;
        
        for i=1:length(syllables);
            Params.SeqFilter.OutlierInds{j}.(syllables{i})=OutlInds{i}; % keep indices
        end
        
        
        
        % REMOVE outliers from data:
        for i=1:length(syllables);
            if ~isempty(OutlInds{i});
                AllDays_RawDatStruct{j}.data.(syllables{i})(OutlInds{i},:)=[]; % remove outlier rows
                disp(['Removed ' num2str(length(OutlInds{i})) ' outliers out of ' num2str(size(AllDays_RawDatStruct{j}.data_WithOutlier.(syllables{i}),1)) ' for syllable ' syllables{i} ' on day ' num2str(j)]);
            end
        end
        
        
        
        
        
        % GET STATS without outliers
        AllDays_RawDatStruct{j}.summary_stats_WithOutlier=AllDays_RawDatStruct{j}.summary_stats; % change name
        
        for i=1:length(syllables);
            X=cell2mat(AllDays_RawDatStruct{j}.data.(syllables{i})(:,1)); % FF vals to matrix
            AllDays_RawDatStruct{j}.summary_stats.(syllables{i}).meanFF=mean(X);
            AllDays_RawDatStruct{j}.summary_stats.(syllables{i}).medianFF=median(X);
            AllDays_RawDatStruct{j}.summary_stats.(syllables{i}).sdFF=std(X);
            AllDays_RawDatStruct{j}.summary_stats.(syllables{i}).n=length(X);
            AllDays_RawDatStruct{j}.summary_stats.(syllables{i}).semFF=std(X)/sqrt(length(X)-1);
            
            % get mean pitch contour
            X=cell2mat(AllDays_RawDatStruct{j}.data.(syllables{i})(:,2)); % rows: rend, cols: timebin
            AllDays_RawDatStruct{j}.summary_stats.(syllables{i}).pitchcountour_mean=mean(X,1);
            
            % get mean spectrogram
            %     X=[];
            %     for ii=1:size(compiled_seqdep_pitch.data_NoOutlier.(syllables{i})(:,4),1); % num rends
            %         X(:,:,ii)=compiled_seqdep_pitch.data_NoOutlier.(syllables{i}){ii,4}; % extract individual specs
            %     end
            %     compiled_seqdep_pitch.summary_stats_NoOutlier.(syllables{i}).spec_mean=mean(X,3);
            
        end
        
        
        % PLOT TO CONFIRM THAT OUTLIERS REMOVED:
        if plotON==1;
            for i = 1:length(syllables);
                if ~isempty(OutlInds{i});
                    
                    %plots the mean and individual renditions
                    figure; hold on, subplot(2,1,1), hold on;
                    PCmat=cell2mat(AllDays_RawDatStruct{j}.data.(syllables{i})(:,2));
                    plot(PCmat','LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
                    plot(mean(PCmat),'Linewidth',2)
                    %         xlabel(['Time, bins: '
                    %         num2str(Params.DayRawDat.pc_T{i}(2)-Params.DayRawDat.pc_T{i}(1))])
                    %         % to include, need to get noteID
                    ylabel('Frequency (Hz)')
                    title(['OUTLIER GONE: Mean and individual pitch contours for ' syllables{i} ' on day ' num2str(j)]);
                    
                    % PLOTS THE OUTLIERS ONLY
                    %plots the mean and individual renditions
                    subplot(2,1,2); hold on;
                    PCmat=cell2mat(AllDays_RawDatStruct{j}.data_WithOutlier.(syllables{i})(OutlInds{i},2));
                    plot(PCmat','LineStyle','--')
                    %         xlabel(['Time, bins: ' num2str(pc_T{i}(2)-pc_T{i}(1))])
                    ylabel('Frequency (Hz)')
                    title(['OUTLIERS THAT WERE REMOVED: ' syllables{i} ' on day ' num2str(j)]);
                end
            end
        end
        
        
        disp(''); disp('');
    end
    
end


%% save params and compiled data

% prepare saving directory
tstamp=lt_get_timestamp(0);
% Params.SeqFilter.savedir=[Params.DayRawDat.savedir '/SeqFilter_' Params.SeqFilter.FirstDay 'to' Params.SeqFilter.LastDay '_' tstamp];
Params.SeqFilter.savedir=[Params.DayRawDat.savedir '/SeqFilterCompile'];

try
    cd(Params.SeqFilter.savedir);
catch err
    
    mkdir(Params.SeqFilter.savedir);
    cd(Params.SeqFilter.savedir);
end

% if old files in dir, moves to OldAnalysis folder
FilesInDir=dir('*.mat');

if length(FilesInDir)>0; 
    
    try 
        cd('OldAnalysis');
        
    catch err
        mkdir('OldAnalysis');
        cd('OldAnalysis');
    end
    
    mkdir(['Moved_' tstamp]);
    cd ../
    
    eval(['!mv * OldAnalysis/Moved_' tstamp]) % moves old save files
end


% save
save('AllDays_RawDatStruct.mat','AllDays_RawDatStruct','-v7.3');
save('Params.mat','Params');


% write a text file that tells you when files were made
fid1=fopen(['DONE_SeqFilterCompile_' tstamp '.txt'],'w');
fclose(fid1);


disp('DONE! seq filtered and outliers removed');





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

%% ---------------------------------------------------------
%% ---------------------------------------------------------
%% TROUBLESHOOTING

%% To plot outliers manually
% % 1) Lists outliers that were found
% AllSyls=fieldnames(AllDays_RawDatStruct{1}.data);
% for i=1:length(AllSyls);
%     
%     syl='b';
%     syl=AllSyls{i};
%     Ndays=size(Params.SeqFilter.OutlierInds,2);
%     % to list outliers detected
%     for i=1:Ndays
%         try
%             OutInds=Params.SeqFilter.OutlierInds{i}.(syl);
%             disp([num2str(length(OutInds)) ' outliers for ' syl ' on day ' num2str(i)]);
%         catch err
%         end
%     end
%     
% end
% 
% % 2) To plot all PCs, overlayed with outliers, time bins, and thresholds
% % (modify threshold params below).
% 
% % enter day and syl and x vals for time bins
% j=18;
% Params.SeqFilter.pc_time_window_list{j}
% Params.DayRawDat.syllables
% 
% x1=42;
% x2=130;
% 
% 
% OutInds=Params.SeqFilter.OutlierInds{j}.(syl);
% 
% % RUN
% figure; hold on;
% % subplot(2,1,1); hold on;
% PCmat=cell2mat(AllDays_RawDatStruct{j}.data.(syl)(:,2));
% plot(PCmat','Color',[0.8, 0.8, 0.8])
% disp(['N(all) = ' num2str(size(PCmat,1))]);
% 
% 
% % subplot(2,1,2); hold on;
% PCmat=cell2mat(AllDays_RawDatStruct{j}.data_WithOutlier.(syl)(OutInds,2));
% 
% plot(PCmat');
% disp(['N(outliers) = ' num2str(size(PCmat,1))]);
% 
% 
% 
% line([x1 x1], ylim);
% line([x2 x2], ylim);
% 
% 
% % plot line showing threshold for outliers
% PCmat=cell2mat(AllDays_RawDatStruct{j}.data_WithOutlier.(syl)(:,2));
% 
% numIQR=2;
% for i=1:size(PCmat,2); % for each column
%     col=PCmat(:,i);
%     high(i)=median(col)+lt_inter90tile(col)/2 + numIQR*lt_inter90tile(col);
%      low(i)=median(col)-lt_inter90tile(col)/2 - numIQR*lt_inter90tile(col);
% end
%     
% % numIQR=1.5;
% % for i=1:size(PCmat,2); % for each column
% %     col=PCmat(:,i);
% %     high(i)=median(col)+iqr(col)/2 + numIQR*iqr(col);
% %      low(i)=median(col)-iqr(col)/2 - numIQR*iqr(col);
% % end
% % 
% 
% plot(1:size(PCmat,2),high,'k');
% plot(1:size(PCmat,2),low,'k');
% 
%     
% 
% % CHECK THAT time windows used for specific syls were correct
% Params.SeqFilter.pc_time_window_list{j}
% AllSyls








