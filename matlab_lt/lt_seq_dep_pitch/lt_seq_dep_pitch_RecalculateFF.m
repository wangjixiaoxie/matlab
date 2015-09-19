function [Params, AllDays_RawDatStruct] = lt_seq_dep_pitch_RecalculateFF(Params, AllDays_RawDatStruct, plotON)
%% LT 5/19/15 - takes new input for time windows for pitch contour, and recalculates FF
% uses data without outliers removed. then removes outliers for the new
% time windows

% Instructions:
% 1) load params and all days stats

% 2) Additional inputs:
% Params.RecalculateFF.pc_time_window_list={[55 100], [20 300]}; % each
% entry is for one syllable. They should be in order of syls in AllDays_RawDatStruct{1}.data

% plotON = 1; plots PCs of all outliers.


%% DISPLAY WARNING

disp('WARNING - this requires the inputed order of time windows to match the order of syllables on day 1!');


%% PARAMS

try
    SylFields=fieldnames(AllDays_RawDatStruct{1}.data);
catch err
    disp('Problem - no data on day 1. I assume syl inds are based on fieldnames of day 1, so cannot continue.')
    keyboard
end

NumSyls= length(SylFields);
NumDays= length(AllDays_RawDatStruct);


%% MAKE SURE ALL DAYS WITH DATA HAVE SAME NUMBER OF SYLS - if not, then some of the windows entered will be incorrect -i.e. they are in order of fieldnames
% note: fieldnames() gets fields in order they were created.  so expected
% to be reliable from day to day.

for i=1:NumDays;
    if ~isempty(AllDays_RawDatStruct{i});
    
    % check that num fields matches num of entries of time windows
    if length(fieldnames(AllDays_RawDatStruct{i}.data)) ~= length(Params.RecalculateFF.pc_time_window_list);
        disp(['Problem, day: ' num2str(i) 'has number of syls not matching other days - code will lead to problem, mismatch window with syl'])
        keyboard
    end
    end
end



%% Recalculate FF
% does this for all syls and using all (outlier-contianing) data. saves the
% data back replacing FF

% go through each day
for ii=1:NumDays;
    if ~isempty(AllDays_RawDatStruct{ii})
        
        % Print sylfields and the actual order of syls. User should
        % verify that they are identical
        disp(['Order of syls for day ' num2str(ii) ' is: ' fieldnames(AllDays_RawDatStruct{ii}.data)'])
        
        % for each syllable
        for i=1:NumSyls;
            syl=SylFields{i};
            
            
            % stats for day
            daystats=AllDays_RawDatStruct{ii}.data_WithOutlier.(syl);
            
            % == pitch contour data
            PCmat=cell2mat(daystats(:,2)); % cell to mat
            
            % what is the PC time window?
            pc_time_window=Params.RecalculateFF.pc_time_window_list(:,i);
            
            % == extract ff
            PCff=mean(PCmat(:,pc_time_window(1):pc_time_window(2)),2); % take mean across time
            
            % == return data to cells
            N=size(daystats,1); % num rends
            AllDays_RawDatStruct{ii}.data_WithOutlier.(syl)(:,1)=mat2cell(PCff,ones(N,1)); % put back into the cells
        end
    end
end


%% Recalculate summary stats (with outliers)
for i=1:NumSyls;
    syl=SylFields{i};
    
    for ii=1:NumDays;
        if ~isempty(AllDays_RawDatStruct{ii});
            
            X=cell2mat(AllDays_RawDatStruct{ii}.data_WithOutlier.(syl)(:,1)); % FF vals to matrix
            
            AllDays_RawDatStruct{ii}.summary_stats_WithOutlier.(syl).meanFF=mean(X);
            AllDays_RawDatStruct{ii}.summary_stats_WithOutlier.(syl).medianFF=median(X);
            AllDays_RawDatStruct{ii}.summary_stats_WithOutlier.(syl).sdFF=std(X);
            AllDays_RawDatStruct{ii}.summary_stats_WithOutlier.(syl).n=length(X);
            AllDays_RawDatStruct{ii}.summary_stats_WithOutlier.(syl).semFF=std(X)/sqrt(length(X)-1);
            
            % get mean pitch contour
            X=cell2mat(AllDays_RawDatStruct{ii}.data_WithOutlier.(syl)(:,2)); % rows: rend, cols: timebin
            AllDays_RawDatStruct{ii}.summary_stats_WithOutlier.(syl).pitchcountour_mean=mean(X,1);
        end
    end
end


%% Refilter outliers, based on the new pitch contour

% == START LOOKING FOR OUTLIERS
for j=1:NumDays;
    
    if ~isempty(AllDays_RawDatStruct{j}); % only continue if day has data.
        
        % -- Start will all data, including outliers (will then throw out outliers)
        AllDays_RawDatStruct{j}.data=AllDays_RawDatStruct{j}.data_WithOutlier;
        
        % -- throw out summary data - they were made with wrong temporal
        % windows.
        if isfield(AllDays_RawDatStruct{j}, 'summary_stats');
            AllDays_RawDatStruct{j}=rmfield(AllDays_RawDatStruct{j},'summary_stats');
        end
        
        % == RUN
        Outliers={};
        OutlInds={};
        for i=1:NumSyls;
            syl=SylFields{i};
            
            % -- make sure this day has this syllable
            if ~isfield(AllDays_RawDatStruct{j}.data, syl);
                continue
            end
            
            % == extract PC for this syl and day
            PC_mat=cell2mat(AllDays_RawDatStruct{j}.data.(syl)(:,2)); % convert to mat
            
            columns_to_check=Params.RecalculateFF.pc_time_window_list(:,i);
            num_columns=columns_to_check(2)-columns_to_check(1)+1; % how many times to do analysis
            Outliers{i}=[];
            
            for ii=1:num_columns; % for each column, perform check.
                % empirically determined 2 to be a good iqr value to use.
                % (based on pu37wh20)
                [~, B, C]=lt_db_tukey_outlier_90tile(PC_mat,columns_to_check(1)+ii-1,2.5);
                %                 [~, B, C]=lt_db_tukey_outlier(PC_mat,columns_to_check(1)+ii-1,1.5);
                
                Outliers{i}=[Outliers{i} B' C']; % accumulate index of outliers
            end
            
            OutlInds{i}=unique(Outliers{i}); % gives back specific index only once
            
            
            % == Annotate which outliers removed for this day and syl
            Params.RecalculateFF.OutlierInds{j}.(syl)=OutlInds{i};
            
            % == Remove outliers from data
            if ~isempty(OutlInds{i});
                AllDays_RawDatStruct{j}.data.(syl)(OutlInds{i},:)=[]; % remove outlier rows
                disp(['Removed ' num2str(length(OutlInds{i})) ' outliers out of ' num2str(size(AllDays_RawDatStruct{j}.data_WithOutlier.(syl),1)) ' for syllable ' syl ' on day ' num2str(j)]);
            end
            
        end
    end
end


% == GET SUMMARY STATS without outliers
for j=1:NumDays;
    
    if isempty(AllDays_RawDatStruct{j})
        continue
    end
    
    for i=1:NumSyls;
        syl=SylFields{i};
        
        if ~isfield(AllDays_RawDatStruct{j}.data, syl); % if this syl has no data today
            continue
        end
        
        X=cell2mat(AllDays_RawDatStruct{j}.data.(syl)(:,1)); % FF vals to matrix
        AllDays_RawDatStruct{j}.summary_stats.(syl).meanFF=mean(X);
        AllDays_RawDatStruct{j}.summary_stats.(syl).medianFF=median(X);
        AllDays_RawDatStruct{j}.summary_stats.(syl).sdFF=std(X);
        AllDays_RawDatStruct{j}.summary_stats.(syl).n=length(X);
        AllDays_RawDatStruct{j}.summary_stats.(syl).semFF=std(X)/sqrt(length(X)-1);
        
        % get mean pitch contour
        X=cell2mat(AllDays_RawDatStruct{j}.data.(syl)(:,2)); % rows: rend, cols: timebin
        AllDays_RawDatStruct{j}.summary_stats.(syl).pitchcountour_mean=mean(X,1);
    end
    
    % PLOT TO CONFIRM THAT OUTLIERS REMOVED:
    if plotON==1;
        for i = 1:NumSyls;
            syl=SylFields{i};
            
            outlinds=Params.RecalculateFF.OutlierInds{j}.(syl);
            
            if ~isempty(outlinds);
                
                %plots the mean and individual renditions
                figure; hold on, subplot(2,1,1), hold on;
                PCmat=cell2mat(AllDays_RawDatStruct{j}.data.(syl)(:,2));
                plot(PCmat','LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
                plot(mean(PCmat),'Linewidth',2)
                
                ylabel('Frequency (Hz)')
                title(['OUTLIER GONE: Mean and individual pitch contours for ' syl ' on day ' num2str(j)]);
                
                % PLOTS THE OUTLIERS ONLY
                %plots the mean and individual renditions
                subplot(2,1,2); hold on;
                PCmat=cell2mat(AllDays_RawDatStruct{j}.data_WithOutlier.(syl)(outlinds,2));
                plot(PCmat','LineStyle','--')
                
                ylabel('Frequency (Hz)')
                title(['OUTLIERS THAT WERE REMOVED: ' syl ' on day ' num2str(j)]);
            end
        end
    end
    
    disp(''); disp('');
end


%% SAVE

% prepare saving directory
tstamp=lt_get_timestamp(0);
Params.RecalculateFF.savedir=[Params.DayRawDat.savedir '/SeqFilterCompile'];

cd(Params.RecalculateFF.savedir);


% save
save('AllDays_RawDatStruct.mat','AllDays_RawDatStruct','-v7.3');
save('Params.mat','Params');

% write a text file that tells you when files were made
fid1=fopen(['DONE_RecalculateFF_' tstamp '.txt'],'w');
fclose(fid1);

disp('DONE! Recalculated FF and removed outliers');




