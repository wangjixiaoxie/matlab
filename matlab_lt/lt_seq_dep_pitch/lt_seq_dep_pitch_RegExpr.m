function [Params, AllDays_RegExpr] = lt_seq_dep_pitch_RegExpr(Params, AllDays_RawDatStruct, saveON, DoLMAN, AllDays_PlotLearning)

%% LT 6/16/15 - take AllDays Raw data structure (i.e. just individual syllables), parse data based on input regular expressions
% Input:
% Params.RegExpr.expressions={'acb+g', 'acb+'};
% AllDays_RawDatStruct, i.e. the output from lt_seq_dep_pitch_SeqFilterCompile
% DoLMAN = 1; then also analyzes LMAN data.  need to have field:
% data_MUSC, e.g. AllDays_RawDatStruct{1}.data_MUSC
% If DoLMAN=1; then need to supply 
% AllDays_PlotLearning (after running lt_seq_dep_pitch_PlotLearning_Musc), since this contains data filtered to within good time window of muscimol activation



%% PARAMS

NumDays=length(AllDays_RawDatStruct);

if ~exist('DoLMAN','var');
    DoLMAN=0;
elseif DoLMAN==1;
    % check that Plot learning struct is entered
    if ~exist('AllDays_PlotLearning' ,'var');
        disp('PROBLEM, need to enter AllDays_PlotLearning as a variable in order to do LMAN analysis)');
        asdceaacvtg; % to halt program.
    end
end

if ~exist('saveON','var');
    saveON=1;
end

Params.PlotRegExpr.AnyNotesWithoutData=0;

%% SINGLE DAYS - For each day in raw dat struct, parse out data - USES DATA THAT HAS NOT HAD OUTLIERS REMOVED

data_field='data_WithOutlier'; % (use just 'data' to use outlier-removed data)
AllDays_RegExpr.day_data=cell(1,NumDays);
for k=1:NumDays;
    
    
    % ======================== PREPARE STUFF
    % Check that today has data
    if isempty(AllDays_RawDatStruct{k})
        continue
    end
    
    % Extract today's data
    RawDatStruct=AllDays_RawDatStruct{k};

    % Prepare output struct for today
    DayData=[];
    
    
    % ============================= RUN
    % First, collect all unique song renditions today and get all syllable
    % sequences/labels
    % - go through all syls and all rends for each syl.  collect song datenums. save the unique ones.
    sylfields=fieldnames(RawDatStruct.(data_field));
    AllSong_datenums=[];
    AllSong_labels={};
    
    for i=1:length(sylfields);
        syl=sylfields{i};
        AllSong_datenums=[AllSong_datenums; cell2mat(RawDatStruct.(data_field).(syl)(:,6))];
        AllSong_labels=[AllSong_labels; RawDatStruct.(data_field).(syl)(:,7)];
    end
    
    [~, IndsToKeep, ~]=unique(AllSong_datenums);
    AllSong_datenums=AllSong_datenums(IndsToKeep); % FINAL, contains unique songs.
    AllSong_labels=AllSong_labels(IndsToKeep);
    
    % ---- OUTPUT DATA
    DayData.AllSong_datenums=AllSong_datenums;
    DayData.AllSong_labels=AllSong_labels;
    
    
    % == Second, Match regular expression - output expected position of single
    % syllables for each detected motif rendition
    reg_expr_list=Params.RegExpr.expressions;
    for j=1:length(reg_expr_list);
        RegExpr=reg_expr_list{j}; % the regular expression to match
        
        % go through all songs and find match
        [StartInds, MatchStrings]=regexp(AllSong_labels, RegExpr, 'start','match');
        
        % ==== IN PROGRESS - through out match syls that are the result of
        % a negation (as they might not have data for that syl.
        
        
        % ============= ORGANIZE OUTPUT ARRAYS AND CELLS
        DayData.(data_field){j}.reg_expr_string=RegExpr;
        DayData.(data_field){j}.StartInds_OfMatch_BySong=StartInds;
        DayData.(data_field){j}.MatchedStrings_BySong=MatchStrings;
        
        
        NumSongs=length(AllSong_datenums);
        motif_counter=1; % all motifs across all songs. (for this regexp)
        
        DayData.(data_field){j}.Final_ARRAYS=[]; % output data goes here
        
        % === COLLECT FF FOR ALL THE MATCHES (by refering back to raw
        % data) (i.e. know song and position for each slot)
        for jj=1:NumSongs;
            song_datenum=AllSong_datenums(jj);
            
            % For this song, go through all motifs.
            % how many motifs for this song?
            NumMotifs=length(StartInds{jj});
            
            for jjj=1:NumMotifs;
                motif_string=MatchStrings{jj}{jjj};
                motif_start_position=StartInds{jj}(jjj);
                
                % for this motif, go through all syls
                for jjjj=1:length(motif_string);
                    
                    % -- Things needed to ID this syllable
                    syl=motif_string(jjjj);
                    syl_pos=motif_start_position+jjjj-1;
                    
                    % =============== TEMPORARY PROBABLY 
                    % if this syllable does not exist in the raw data, put
                    % a nan in its place
                    if ~isfield(RawDatStruct.(data_field), syl);
                        % then does not exist
                        
                        DayData.(data_field){j}.Final_ARRAYS.FFvals(motif_counter,jjjj)=nan;
                        DayData.(data_field){j}.Final_ARRAYS.Trig(motif_counter,jjjj)=nan;
                        DayData.(data_field){j}.Final_ARRAYS.CatchTrial(motif_counter,jjjj)=nan;
                        
                        DayData.(data_field){j}.Final_ARRAYS.RawDat_SylPos(motif_counter,jjjj)=syl_pos;
                        DayData.(data_field){j}.Final_ARRAYS.RawDatInd_WithOutliers(motif_counter,jjjj)=nan;
                        
                        % Make a note in params
                        Params.PlotRegExpr.AnyNotesWithoutData=1;
                        
                        continue;
                    end
                    
                    % -- Gather information about this syl
                    potential_inds_bysong=find(cell2mat(RawDatStruct.(data_field).(syl)(:,6))==song_datenum); % find this song
                    potential_inds_bypos=find(cell2mat(RawDatStruct.(data_field).(syl)(:,8))==syl_pos); % find this position
                    
                    Ind_of_data=intersect(potential_inds_bysong, potential_inds_bypos);
                    
                    % ASIDE ---- make sure the detected ind is length 1
                    if isempty(Ind_of_data)
                        disp('WARNING - a case where found no data to slot into reg_expr detected data');
                    elseif length(Ind_of_data)>1;
                        disp('WARNING - a case where found >1 datapoint to slot into reg_expr detected data');
                    end
                    % ------------------
                    
                    % === OUTPUT DATA PER SYLLABLE
                    FFval=RawDatStruct.(data_field).(syl){Ind_of_data,1};
                    Trig=RawDatStruct.(data_field).(syl){Ind_of_data,9};
                    CatchTrial=RawDatStruct.(data_field).(syl){Ind_of_data,9};
                    
                    DayData.(data_field){j}.Final_ARRAYS.FFvals(motif_counter,jjjj)=FFval;
                    DayData.(data_field){j}.Final_ARRAYS.Trig(motif_counter,jjjj)=Trig;
                    DayData.(data_field){j}.Final_ARRAYS.CatchTrial(motif_counter,jjjj)=CatchTrial;
                    
                    DayData.(data_field){j}.Final_ARRAYS.RawDat_SylPos(motif_counter,jjjj)=syl_pos;
                    DayData.(data_field){j}.Final_ARRAYS.RawDatInd_WithOutliers(motif_counter,jjjj)=Ind_of_data;
                    
                end
                
                % === OUTPUT DATA FOR THIS MOTIF
                DayData.(data_field){j}.Final_ARRAYS.SylStrings{motif_counter}=motif_string;
                DayData.(data_field){j}.Final_ARRAYS.Song_DateNum(motif_counter)=song_datenum;
                
                % --- iterate one motif
                motif_counter=motif_counter+1;
            end
        end
        
        % Reshape for fun
        DayData.(data_field){j}.Final_ARRAYS.SylStrings=DayData.(data_field){j}.Final_ARRAYS.SylStrings';
        DayData.(data_field){j}.Final_ARRAYS.Song_DateNum=DayData.(data_field){j}.Final_ARRAYS.Song_DateNum';
 

        % === Find empty positions in array and replace with nan
        % use fact that a real FF cannot be 0
        Inds_To_Remove=DayData.(data_field){j}.Final_ARRAYS.FFvals==0;
        
        DayData.(data_field){j}.Final_ARRAYS.FFvals(Inds_To_Remove)=nan;
        DayData.(data_field){j}.Final_ARRAYS.Trig(Inds_To_Remove)=nan;
        DayData.(data_field){j}.Final_ARRAYS.CatchTrial(Inds_To_Remove)=nan;
        DayData.(data_field){j}.Final_ARRAYS.RawDat_SylPos(Inds_To_Remove)=nan;
        DayData.(data_field){j}.Final_ARRAYS.RawDatInd_WithOutliers(Inds_To_Remove)=nan;
              
    end
    
    % ======= SLIDE THIS DAYS DATA INTO ALL DAYS STRUCTURE
    AllDays_RegExpr.day_data{k}=DayData;
end

%% [REPEAT ANALYSIS FOR MUSCIMOL] SINGLE DAYS - For each day in raw dat struct, parse out data - USES DATA THAT HAS NOT HAD OUTLIERS REMOVED

if DoLMAN==1;
    
    data_field='data_WithOutlier';
    
    syls_unique=Params.PlotLearning.SylFields_Unique;
    
    % ================== FIGURE OUT GOOD DAYS AND GOOD SONGS
    MuscimolDayInformation.MuscimolDays_good=find(~cellfun(@isempty, AllDays_PlotLearning.DataMatrix_MUSC.(syls_unique{1}).Tvals_WithinTimeWindow));
    
    % for each day, get the Tvals of the good songs
    for i=1:length(MuscimolDayInformation.MuscimolDays_good);
        day=MuscimolDayInformation.MuscimolDays_good(i);
        
        % -- Go through all syls and collect song information and song label
        % information
        tvals_within_time_window_ALL=[];
        song_labels_ALL={};
        
        for ii=1:length(syls_unique);
            syl=syls_unique{ii};
            
            tvals_thissyl=AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
            tvals_thissyl=unique(tvals_thissyl);
            
            % with those tvals, go to raw dat and get song labels
            labels_thissyl={};
            for iii=1:length(tvals_thissyl);
                tval=tvals_thissyl(iii);
                
                % for this tval, find the corresponging song label
                song_inds=find(cell2mat(AllDays_RawDatStruct{day}.data_MUSC.(syl)(:,6))==tval);
                
                % --- ASIDE: If don't find the song ind, then there is problem
                if isempty(song_inds);
                    disp('Note: tval should be oen of days songs but cannot find it..., assuming is because took MUSC data outside of musc window');
                    
                    %  keyboard;
                    % then is likely because this datapoint was taken from
                    % PBS, in the window after musc.
                    song_inds=find(cell2mat(AllDays_RawDatStruct{day}.data.(syl)(:,6))==tval);
                    labels_thissyl{iii}=AllDays_RawDatStruct{day}.data.(syl){song_inds(1),7};
                    
                else
                    % pick out the label for one of the song inds
                    labels_thissyl{iii}=AllDays_RawDatStruct{day}.data_MUSC.(syl){song_inds(1),7};
                end
            end
            
            if isempty(labels_thissyl)
                disp('ACTUAL PROBLEM - cant find this song using tvals');
            end
            
            
            
            % ---- Put tvals and labels into output variables
            tvals_within_time_window_ALL=[tvals_within_time_window_ALL tvals_thissyl];
            song_labels_ALL=[song_labels_ALL labels_thissyl];
        end
        
        % GET UNIQUE SONGS (and corresponding song labels)
        [tvals_within_time_window_ALL, inds, ~]=unique(tvals_within_time_window_ALL);
        song_labels_ALL=song_labels_ALL(inds);
        
        % === OUTPUT INFO
        MuscimolDayInformation.AllSongs_WithinTimeWindow{day}=tvals_within_time_window_ALL;
        MuscimolDayInformation.AllLabels_WithinTimeWindow{day}=song_labels_ALL;
    end
    
    
    % =================== PERFORM ANALYSIS
    AllDays_RegExpr.day_data_MUSC=cell(1,NumDays);
    for n=1:length(MuscimolDayInformation.MuscimolDays_good);
        
        k=MuscimolDayInformation.MuscimolDays_good(n); % day ind
        
        % ======================== PREPARE STUFF
        % Check that today has data
        if isempty(AllDays_RawDatStruct{k})
            continue
        end
        
        % Extract today's data
        RawDatStruct=AllDays_RawDatStruct{k};
        
        % Prepare output struct for today
        DayData=[];
        
        
        % ============================= RUN
        % ========== First, collect all unique song renditions today and get all syllable
        % sequences/labels
        % 1) Datenums are the ones extracted for in good muscimol window
        AllSong_datenums=unique(MuscimolDayInformation.AllSongs_WithinTimeWindow{k});
        
        % 2) For those datenums, go find the song labels
        AllSong_labels=MuscimolDayInformation.AllLabels_WithinTimeWindow{k};
        
        % ---- OUTPUT DATA
        DayData.AllSong_datenums=AllSong_datenums;
        DayData.AllSong_labels=AllSong_labels;
        % =============================
        
        % == Second, Match regular expression - output expected position of single
        % syllables for each detected motif rendition
        reg_expr_list=Params.RegExpr.expressions;
        for j=1:length(reg_expr_list);
            RegExpr=reg_expr_list{j}; % the regular expression to match
            
            % go through all songs and find match
            [StartInds, MatchStrings]=regexp(AllSong_labels, RegExpr, 'start','match');
            
            
            % ==== IN PROGRESS - through out match syls that are the result of
            % a negation (as they might not have data for that syl.
            
            
            % ============= ORGANIZE OUTPUT ARRAYS AND CELLS
            DayData.(data_field){j}.reg_expr_string=RegExpr;
            DayData.(data_field){j}.StartInds_OfMatch_BySong=StartInds;
            DayData.(data_field){j}.MatchedStrings_BySong=MatchStrings;
            
            % skip if don't match at all (no song has data)
            tmp=[];
            for ll=1:length(StartInds);
                if ~isempty(StartInds{ll});
                tmp=1;
                end
            end
            if isempty(tmp)
                continue;
            end
            
            % Continue
            NumSongs=length(AllSong_datenums);
            motif_counter=1; % all motifs across all songs. (for this regexp)
            
            DayData.(data_field){j}.Final_ARRAYS=[]; % output data goes here
            
            % === COLLECT FF FOR ALL THE MATCHES (by refering back to raw
            % data) (i.e. know song and position for each slot)
            for jj=1:NumSongs;
                song_datenum=AllSong_datenums(jj);
                
                % For this song, go through all motifs.
                % how many motifs for this song?
                NumMotifs=length(StartInds{jj});
                
                if NumMotifs>0;
                    for jjj=1:NumMotifs;
                        motif_string=MatchStrings{jj}{jjj};
                        motif_start_position=StartInds{jj}(jjj);
                        
                        % for this motif, go through all syls
                        for jjjj=1:length(motif_string);
                            
                            % -- Things needed to ID this syllable
                            syl=motif_string(jjjj);
                            syl_pos=motif_start_position+jjjj-1;
                            
                            
                            % =============== TEMPORARY PROBABLY
                            % if this syllable does not exist in the raw data, put
                            % a nan in its place
                            if ~isfield(RawDatStruct.data_WithOutlier, syl);
                                % then does not exist
                                
                                DayData.(data_field){j}.Final_ARRAYS.FFvals(motif_counter,jjjj)=nan;
                                DayData.(data_field){j}.Final_ARRAYS.Trig(motif_counter,jjjj)=nan;
                                DayData.(data_field){j}.Final_ARRAYS.CatchTrial(motif_counter,jjjj)=nan;
                                
                                DayData.(data_field){j}.Final_ARRAYS.RawDat_SylPos(motif_counter,jjjj)=syl_pos;
                                DayData.(data_field){j}.Final_ARRAYS.RawDatInd_WithOutliers(motif_counter,jjjj)=nan;
                                
                                Params.PlotRegExpr.AnyNotesWithoutData=1;
                                continue;
                            end
                            
                            % -- Gather information about this syl
                            potential_inds_bysong=find(cell2mat(RawDatStruct.data_WithOutlier.(syl)(:,6))==song_datenum); % find this song
                            potential_inds_bypos=find(cell2mat(RawDatStruct.data_WithOutlier.(syl)(:,8))==syl_pos); % find this position
                            
                            Ind_of_data=intersect(potential_inds_bysong, potential_inds_bypos);
                            
                            % ASIDE ---- make sure the detected ind is length 1
                            if isempty(Ind_of_data)
                                disp('WARNING - a case where found no data to slot into reg_expr detected data');
                            elseif length(Ind_of_data)>1;
                                disp('WARNING - a case where found >1 datapoint to slot into reg_expr detected data');
                            end
                            % ------------------
                            
                            % === OUTPUT DATA PER SYLLABLE
                            FFval=RawDatStruct.data_WithOutlier.(syl){Ind_of_data,1};
                            Trig=RawDatStruct.data_WithOutlier.(syl){Ind_of_data,9};
                            CatchTrial=RawDatStruct.data_WithOutlier.(syl){Ind_of_data,9};
                            
                            DayData.(data_field){j}.Final_ARRAYS.FFvals(motif_counter,jjjj)=FFval;
                            DayData.(data_field){j}.Final_ARRAYS.Trig(motif_counter,jjjj)=Trig;
                            DayData.(data_field){j}.Final_ARRAYS.CatchTrial(motif_counter,jjjj)=CatchTrial;
                            
                            DayData.(data_field){j}.Final_ARRAYS.RawDat_SylPos(motif_counter,jjjj)=syl_pos;
                            DayData.(data_field){j}.Final_ARRAYS.RawDatInd_WithOutliers(motif_counter,jjjj)=Ind_of_data;
                            
                        end
                        
                        % === OUTPUT DATA FOR THIS MOTIF
                        DayData.(data_field){j}.Final_ARRAYS.SylStrings{motif_counter}=motif_string;
                        DayData.(data_field){j}.Final_ARRAYS.Song_DateNum(motif_counter)=song_datenum;
                        
                        % --- iterate one motif
                        motif_counter=motif_counter+1;
                    end
                end
            end
            
            % Reshape for fun
            DayData.(data_field){j}.Final_ARRAYS.SylStrings=DayData.(data_field){j}.Final_ARRAYS.SylStrings';
            DayData.(data_field){j}.Final_ARRAYS.Song_DateNum=DayData.(data_field){j}.Final_ARRAYS.Song_DateNum';
            
            % convert datenum to tvals (hours)
            
            
            % === Find empty positions in array and replace with nan
            % use fact that a real FF cannot be 0
            Inds_To_Remove=DayData.(data_field){j}.Final_ARRAYS.FFvals==0;
            
            DayData.(data_field){j}.Final_ARRAYS.FFvals(Inds_To_Remove)=nan;
            DayData.(data_field){j}.Final_ARRAYS.Trig(Inds_To_Remove)=nan;
            DayData.(data_field){j}.Final_ARRAYS.CatchTrial(Inds_To_Remove)=nan;
            DayData.(data_field){j}.Final_ARRAYS.RawDat_SylPos(Inds_To_Remove)=nan;
            DayData.(data_field){j}.Final_ARRAYS.RawDatInd_WithOutliers(Inds_To_Remove)=nan;
            
        end
        
        % ======= SLIDE THIS DAYS DATA INTO ALL DAYS STRUCTURE
        AllDays_RegExpr.day_data_MUSC{k}=DayData;
    end
    
    Params.RegExpr.MuscimolDayInformation=MuscimolDayInformation;
    
end

%% Calculate motif probabilities across days

for i=1:NumDays;
    
end


%% PARSE DATA based on the actual subclass of the regular expression
% e.g. acb+d could be acbd, acbbd, acbbbbd, etc.  
data_field='data_WithOutlier';

% ===== GO THROUGH all reg expr strings, and for each string,
% collect all the possible subclasses across all days
Params.RegExpr.subexpressions={};

for i=1:length(Params.RegExpr.expressions);
    regexpr_string=Params.RegExpr.expressions{i};
     
    % ============================ Across all days, get list of all classes that exist in at
    % least one day
    
    AllSubClasses=[];
    for j=1:NumDays;
        if  isempty(AllDays_RegExpr.day_data{j});
            continue
        end
        AllSubClasses=[AllSubClasses; AllDays_RegExpr.day_data{j}.(data_field){i}.Final_ARRAYS.SylStrings];
    end
    
    AllSubClasses=unique(AllSubClasses);
    
    % --- Sort by order of length
    tmp_lengths=cellfun(@length,AllSubClasses);
    [~, inds]=sort(tmp_lengths);
    
    AllSubClasses=AllSubClasses(inds);
    
    % --- Put into params
    Params.RegExpr.subexpressions{i}=AllSubClasses;
    
    % ============================= SLOT INTO OUTPUT STRUCTURE
    for j=1:length(AllSubClasses);
        sub_class=AllSubClasses{j};
        
        for jj=1:NumDays;
            
            if isempty(AllDays_RegExpr.day_data{jj});
                continue
            end
            
            % -- find renditions this day that are in this repeat class
            AllSylStrings=AllDays_RegExpr.day_data{jj}.(data_field){i}.Final_ARRAYS.SylStrings;
            
            Inds_ThisRepClass=strcmp(AllSylStrings, sub_class);
            
            % ========= GATHER DAY DATA FOR THIS SUBCLASS
            FFvals_All=AllDays_RegExpr.day_data{jj}.(data_field){i}.Final_ARRAYS.FFvals;
            FFvals_ThisRepClass=FFvals_All(Inds_ThisRepClass,:);
            
            Tvals_All=AllDays_RegExpr.day_data{jj}.(data_field){i}.Final_ARRAYS.Song_DateNum;
            Tvals_ThisRepClass=Tvals_All(Inds_ThisRepClass);
            
            SylPos_All=AllDays_RegExpr.day_data{jj}.(data_field){i}.Final_ARRAYS.RawDat_SylPos;
            SylPos_ThisRepClass=SylPos_All(Inds_ThisRepClass);
            
            % === OUTPUT DATA
            AllDays_RegExpr.day_data{jj}.data_ParsedIntoSubclasses{i}.regexpr_string=regexpr_string;
            AllDays_RegExpr.day_data{jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.string=sub_class;
            
            AllDays_RegExpr.day_data{jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.FFvals=FFvals_ThisRepClass;
            
            AllDays_RegExpr.day_data{jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.Tvals=Tvals_ThisRepClass;
            AllDays_RegExpr.day_data{jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.SylPos=SylPos_ThisRepClass;
            AllDays_RegExpr.day_data{jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.Inds_InUnparsedData=Inds_ThisRepClass;
            
        end
    end
end


%% [MUSCIMOL] PARSE DATA based on the actual subclass of the regular expression
% e.g. acb+d could be acbd, acbbd, acbbbbd, etc.

if DoLMAN==1;
    data_field='data_WithOutlier';
    day_data_field='day_data_MUSC';
    
    % ===== GO THROUGH all reg expr strings, and for each string,
    % collect all the possible subclasses across all days    
    for i=1:length(Params.RegExpr.expressions);
        regexpr_string=Params.RegExpr.expressions{i};
        
        % ============================ Across all days, get list of all classes that exist in at
        % least one day
        AllSubClasses=Params.RegExpr.subexpressions{i};
        
        % ============================= SLOT INTO OUTPUT STRUCTURE
        for j=1:length(AllSubClasses);
            sub_class=AllSubClasses{j};
            
            for jj=MuscimolDayInformation.MuscimolDays_good;
                                
                if isempty(AllDays_RegExpr.(day_data_field){jj});
                    continue
                end
                
                if ~isfield(AllDays_RegExpr.(day_data_field){jj}.(data_field){i}, 'Final_ARRAYS');
                    continue
                end
                
                % -- find renditions this day that are in this repeat class
                AllSylStrings=AllDays_RegExpr.(day_data_field){jj}.(data_field){i}.Final_ARRAYS.SylStrings;
                
                Inds_ThisRepClass=strcmp(AllSylStrings, sub_class);
                
                % ========= GATHER DAY DATA FOR THIS SUBCLASS
                FFvals_All=AllDays_RegExpr.(day_data_field){jj}.(data_field){i}.Final_ARRAYS.FFvals;
                FFvals_ThisRepClass=FFvals_All(Inds_ThisRepClass,:);
                
                Tvals_All=AllDays_RegExpr.(day_data_field){jj}.(data_field){i}.Final_ARRAYS.Song_DateNum;
                Tvals_ThisRepClass=Tvals_All(Inds_ThisRepClass);
                
                SylPos_All=AllDays_RegExpr.(day_data_field){jj}.(data_field){i}.Final_ARRAYS.RawDat_SylPos;
                SylPos_ThisRepClass=SylPos_All(Inds_ThisRepClass);
                
                % === OUTPUT DATA
                AllDays_RegExpr.(day_data_field){jj}.data_ParsedIntoSubclasses{i}.regexpr_string=regexpr_string;
                AllDays_RegExpr.(day_data_field){jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.string=sub_class;
                
                AllDays_RegExpr.(day_data_field){jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.FFvals=FFvals_ThisRepClass;
                
                AllDays_RegExpr.(day_data_field){jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.Tvals=Tvals_ThisRepClass;
                AllDays_RegExpr.(day_data_field){jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.SylPos=SylPos_ThisRepClass;
                AllDays_RegExpr.(day_data_field){jj}.data_ParsedIntoSubclasses{i}.sub_class{j}.Inds_InUnparsedData=Inds_ThisRepClass;
                
            end
        end
    end
    
end

%% SAVE
if saveON==1;
    timestampSv=lt_get_timestamp(0);
    cd(Params.SeqFilter.savedir);
    
    save('Params','Params');
    save('AllDays_RegExpr','AllDays_RegExpr');
    
    % write a text file that tells you when files were made
    fid1=fopen(['DONE_RegExpr_' timestampSv '.txt'],'w');
    fclose(fid1);

end



