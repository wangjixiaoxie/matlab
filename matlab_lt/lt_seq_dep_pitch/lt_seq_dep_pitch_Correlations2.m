function lt_seq_dep_pitch_Correlations(Params, AllDays_RawDatStruct, DaysWanted)
%% LT 6/9/15 - new version, branching from lt_seq_dep_pitch_Correlations;
% Here will use raw data and parse songs into things like motifs, etc.

% INSTRUCTIONS:
% First run Plot Learning.

%% NOTES BELOW: from lt_seq_dep_pitch_Correlations

%% LT 4/11/15 - Look at correlations between syllables
% TO DO:
% 1) How does correlation change thru learnig? - shoudl regress out
% learning. does that change relate to change in generalization?
% 2) Plot learning vs. correaltions (and structural similarity) - 
% 3) plot scatter of trials to get sense of strength of correlation.
% 4) Make sure corr and/or z-scoring is correct. got corr of 100 for using
% z-scores.
% 5) check that using song is a good proxy
% 6) plot sample sizes, etc.
% 7) seems to throw out a lot of songs. check that they are actually
% missing that data

% Inputs
% DaysWanted = 'baseline' % gets baseline corrs
% DaysWanted =10:13 % days 10:13


%% PARAMS

% syllables
SylFields_Unique=Params.PlotLearning.SylFields_Unique;


%% PLOT ALL RENDITIONS RAW DATA TO EYEBALL CORRELATIONS
% one set of plots for each baseline day

for j=1:length(Params.SeqFilter.BaselineDays);
    day=Params.SeqFilter.BaselineDays(j);
    if ~isempty(AllDays_RawDatStruct{day})
        for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder); % for each motif
            
            lt_figure; hold on;
            syllist=Params.SeqFilter.SylLists.FieldsInOrder{i};
            
            for ii=1:length(syllist); % for each syllable
                syl=syllist{ii};
                
                lt_subplot(length(syllist), 1, ii); hold on;
                ylabel(syl);
                
                % === PLOT all rends for this syl on this day
                TT=cell2mat(AllDays_RawDatStruct{day}.data.(syl)(:,6));
                [~, TT]=lt_convert_datenum_to_hour(TT);
                Tvals=TT.hours;
                
                FFvals=cell2mat(AllDays_RawDatStruct{day}.data.(syl)(:,1));
                
                % Plot deviation from day mean as bar plot
                FFvals_dev_day_mean=FFvals-mean(FFvals);
                
                % change color depending on deviation sign
                Inds_positive=FFvals_dev_day_mean>0;
                Inds_negative=FFvals_dev_day_mean<0;
                
                % ----PLOT
                %             stem(Tvals(Inds_positive), FFvals_dev_day_mean(Inds_positive),'Color','b');
                %             stem(Tvals(Inds_negative), FFvals_dev_day_mean(Inds_negative),'Color','r');
                
                lt_plot(Tvals(Inds_positive), FFvals_dev_day_mean(Inds_positive),{'Color','b'});
                lt_plot(Tvals(Inds_negative), FFvals_dev_day_mean(Inds_negative),{'Color','r'});
                
                lt_plot_zeroline;
                
                % remove x label
                if ii~=length(syllist);
                    set(gca, 'XTick',[]);
                end
            end
            
            lt_subtitle(['day: ' num2str(day)]);

        end
    end
end
    
    



%% CORRELATIONS BETWEE SYLLABLES IN THE SAME MOTIF - trial to trial
% Every motif will have a matrix - columns are syls in order, and rows are
% all trials

%% WARNINGS TO USER:
disp('WARNING - correlations currently only works if motifs in Params are actual motifs (i.e. no gaps)')

% What are the motifs?
MotifList=Params.SeqFilter.SylLists.FieldsInOrder;

for i=1:length(MotifList);
    motif_syls=MotifList{i};
    
    
    
    % == Make an array for this motif for this day. also have a
    % corresponding vectors for 1) song number, and 2) position in
    % song (array, one col for each syl); one row per motif
    % STRATEGY - choose the first syl in motif. since all motif rends
    % should have this syl, go through all rends of this syl. for each
    % one store FF and record the song number and the position in song.
    % then go through all the syls in order and slot values into the
    % corresponding spot. keep motifs even if they are incomplete.
    % === SEED starting arrays using the first syl in each motif.
    %         NumRends=size(AllDays_RawDatStruct{day}.data.(syl),1);
    
    first_syl=motif_syls{1};
    
    motif_array_AllDays=[];
    song_num_array_AllDays=[];
    motifbegin_pos_in_song_AllDays=[];
    
    
    for ii=1:length(Params.SeqFilter.BaselineDays);
        day=Params.SeqFilter.BaselineDays(ii);
        
        if isempty(AllDays_RawDatStruct{day});
            continue
        end
        
        % --- Output arrays for this day.
        motif_array=[];
        song_num_array=[];
        motifbegin_pos_in_song=[];
        
        % --- seed the arrays
        motif_array(:,1)=cell2mat(AllDays_RawDatStruct{day}.data.(first_syl)(:,1)); % FF vals
        song_num_array=cell2mat(AllDays_RawDatStruct{day}.data.(first_syl)(:,6)); % datenums
        motifbegin_pos_in_song(:,1)=cell2mat(AllDays_RawDatStruct{day}.data.(first_syl)(:,8)); % note position
        
        
        % === Go through all other syls in this motif and stick data into
        % the correct location
        for j=2:length(motif_syls); % start from 2 since 1 is first syl (above)
            syl=motif_syls{j};
            
            % for each rendition find correct slot in arrays and stick in data
            NumRends=size(AllDays_RawDatStruct{day}.data.(syl),1);
            
            
            for jj=1:NumRends;
                
                % --- GET DATA for this syl
                FFval=AllDays_RawDatStruct{day}.data.(syl){jj,1};
                songnum=AllDays_RawDatStruct{day}.data.(syl){jj,6};
                posinsong=AllDays_RawDatStruct{day}.data.(syl){jj,8};
                
                % ===== FIND WHERE TO SLOT THESE DATA
                % based on song number and motif number,
                % what position in the array should we put this data?
                row_to_slot_data=[];
                
                potential_rows_to_slot_data=find(song_num_array==songnum); % will be as many as there are motifs in this song
                
                % rarely this song will not exist - i.e. motif missing
                % first syl in a song. in that case start new row
                if isempty(potential_rows_to_slot_data);
                    % then put a new row with this data.
                    row_to_slot_data=length(song_num_array)+1;
                    
                else
                    % Found this datapoint's song, which row does it go into?
                    % a song can be mutliple rows. so which row is the
                    % correct motif?
                    
                    % for each motif that is found for this song (i.e for each row), check
                    % to make sure the position of the first syllable is what one would
                    % expect
                    
                    positions_of_motif_starts=motifbegin_pos_in_song(potential_rows_to_slot_data,1);
                    
                    % -- this syl should be a specific number of notes away from that
                    % start
                    ind_temp=positions_of_motif_starts==posinsong-j+1;
                    row_to_slot_data=potential_rows_to_slot_data(ind_temp);
                    
                    % ----- Potential problem, none of motifs is
                    % correct. in that case could be because the first syl is missing, but other syls are there
                    if isempty(row_to_slot_data);
                        all_syl_start_positions=motifbegin_pos_in_song(potential_rows_to_slot_data,:); % for all potential rows, get all syls that have already been slotted in (i.e. 2d array)
                        if j>2
                            for k=2:j-1; % go through all syls preceding the current one (skipping first, as already did)
                                ind_temp=all_syl_start_positions==posinsong-j+k;
                                ind_temp=sum(ind_temp,2);
                                
                                if sum(ind_temp)>0;
                                    % then we have a row to slot into
                                    ind_temp=logical(ind_temp);
                                    row_to_slot_data=potential_rows_to_slot_data(ind_temp);
                                    break
                                end
                            end
                        end
                        
                        if isempty(row_to_slot_data); % if is STILL empty
                            
                            % motif does not already have a row
                            % start a new row.
                            row_to_slot_data=length(song_num_array)+1;
                        end
                        
                    end
                    
                end
                
                % --- Problem - more than one row
                if length(row_to_slot_data)>1;
                    disp('PROBLEM - more than one row found...');
                end
                
                % === FINALLY, SLOT in data
                motif_array(row_to_slot_data,j)=FFval;
                song_num_array(row_to_slot_data)=songnum; % datenums
                motifbegin_pos_in_song(row_to_slot_data,j)=posinsong; % note position
            end
        end
        
        
        % === make sure arrays are correct size and match each other.
        if size(motif_array,1)~=length(song_num_array) || size(motif_array,1)~=size(motifbegin_pos_in_song,1);
            disp('PROBLEM WITH FINAL MOTIF SIZES');
        end
        
        if size(motif_array,2)~=size(motifbegin_pos_in_song,2) || size(motif_array,2)~=length(motif_syls);
            disp('WRONG NUMBER OF SYLS (e.g. a, b, ...) IN ARRAY')
        end
        
        % === SORT ARRAYS by date and motif num
        [song_num_array, inds]=sort(song_num_array);
        motif_array=motif_array(inds,:);
        motifbegin_pos_in_song=motifbegin_pos_in_song(inds,:);
        
        % === REPLACE 0 vals with nan.
        % note that any 0 cannot be real data (can't be FF, note position, or
        % datenum) - represents lack of data.
        motif_array(motif_array==0)=nan;
        motifbegin_pos_in_song(motifbegin_pos_in_song==0)=nan;
        
        
        % ======================= ADD THIS DAYS DATA TO PREVIOUS DAYS
        motif_array_AllDays=[motif_array_AllDays; motif_array];
        song_num_array_AllDays=[song_num_array_AllDays; song_num_array];
        motifbegin_pos_in_song_AllDays=[motifbegin_pos_in_song_AllDays; motifbegin_pos_in_song];
        
    end
    
    
    % =========== STORE MOTIF DATA (ACROSS ALL DAYS) IN OUTPUT STRUCTURE
    AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.Song_datenum_vals=song_num_array_AllDays;
    AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.FFvals_array_containsNaN=motif_array_AllDays;
    AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.SylPos_in_song=motifbegin_pos_in_song_AllDays;
    AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.SylsInMotif=motif_syls;
end



% === VERIFY THINGS ARE CORRECT
disp(' ');
for i=1:length(MotifList);
    
    motif_syls=MotifList{i};
    
    if any(any(diff(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.SylPos_in_song,1,2)>1));
        disp('problem: found syls not in order');
    end
    
    % how many songs?
    unique_songs=unique(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.Song_datenum_vals(:));
    disp([num2str(length(unique_songs)) ' songs for motif num: ' num2str(i)]);
    
    num_motifs=length(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.Song_datenum_vals(:));
    disp([num2str(num_motifs) ' motifs']);
    
    % check that mean pitch and rends and all that match what expected
    
    for ii=1:length(motif_syls);
        syl=motif_syls{ii};
        
        % === get stats from the motif sorted data above
        FFvals=AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.FFvals_array_containsNaN(:,ii);
        songnums=AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.Song_datenum_vals;
        sylposvals=AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.SylPos_in_song(:,ii);
        
        
        inds_to_throwout=isnan(FFvals);
        FFvals(inds_to_throwout)=[];
        songnums(inds_to_throwout)=[];
        sylposvals(inds_to_throwout)=[];
        
        
        N=length(FFvals);
        FFmean=mean(FFvals);
        FFstd=std(FFvals);
        
        disp(['Syllable: ' syl]);
        disp(['MOTIF CALCULATED STATS:'])
        disp(['N=' num2str(N) '; FFmean=' num2str(FFmean) '; FFstd=' num2str(FFstd)]) 
       
        
        % === get previously calcualted stats and print for user
        N_prev=AllDays_PlotLearning.EpochData.Baseline.(syl).n;
        FFmean_prev=AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
        FFstd_prev=AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
        
        disp(['PREVIOUSLY CALCULATED STATS:'])
        disp(['N=' num2str(N_prev) '; FFmean=' num2str(FFmean_prev) '; FFstd=' num2str(FFstd_prev)]) 
        disp(' ');
        
        % ==== PLOT across time to verify that it is correct
        lt_figure; hold on;
        
        % ---- MOTIF STATS
        
        h(1)=lt_subplot(3,1,1); hold on; title(syl);
        plot(FFvals,'or');
        
        h(2)=lt_subplot(3,1,2); hold on;
        plot(sylposvals,'or');
        
        h(3)=lt_subplot(3,1,3); hold on;
        plot(songnums,'or');
        
        linkaxes(h,'x')
        pause;
        
        % === PLOT RAW DATA IN SAME WAY TO COMPARE
        FFvals=[];
        songnums=[];
        sylposvals=[];
        for day=Params.SeqFilter.BaselineDays
            FFvals=[FFvals; cell2mat(AllDays_RawDatStruct{day}.data.(syl)(:,1))];
            songnums=[songnums; cell2mat(AllDays_RawDatStruct{day}.data.(syl)(:,6))];
            sylposvals=[sylposvals; cell2mat(AllDays_RawDatStruct{day}.data.(syl)(:,8))];
        end
        
        % -- PLOT
        
        h(1)=lt_subplot(3,1,1); hold on; title(syl);
        plot(1:length(FFvals), FFvals,'ob');
        
        h(2)=lt_subplot(3,1,2); hold on;
        plot(1:length(sylposvals), sylposvals,'ob');
        
        h(3)=lt_subplot(3,1,3); hold on;
        plot(1:length(songnums), songnums,'ob');
        
        linkaxes(h,'x')
        
        % ----------
        pause
        close all;

       
        
    end
    
end

% ===== SOME OUTPUT PARAMS
NumMotifs=length(MotifList);


%% TROUBLESHOOTING

% motif_num=1;
% syl='aB';
% syl_pos=2;
% 
% % go through songs one by one. print positions and FFvals
% unique_songs=unique(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{motif_num}.Song_datenum_vals(:));
%     
% rawsongs_inorder=cell2mat(AllDays_RawDatStruct{1}.data.(syl)(:,6));
% 
% for i=1:length(unique_songs);
%     % motif_num=1;
% syl='aB';
% syl_pos=2;
% 
% % go through songs one by one. print positions and FFvals
% unique_songs=unique(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{motif_num}.Song_datenum_vals(:));
%     
% rawsongs_inorder=cell2mat(AllDays_RawDatStruct{1}.data.(syl)(:,6));
% 
% for i=1:length(unique_songs);
%     
%     songnum=unique_songs(i);
%     
%     find(rawsongs_inorder==songnum);
%     
% end
% 
% figure; hold on;
% plot(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{motif_num}.FFvals_array_containsNaN(:,syl_pos));
% 
% % == raw data
% 
% rawdates=cell2mat(AllDays_RawDatStruct{1}.data.(syl)(:,
% plot(cell2mat(AllDays_RawDatStruct{1}.data.(syl)(:,1)))

%     songnum=unique_songs(i);
%     
%     find(rawsongs_inorder==songnum);
%     
% end
% 
% figure; hold on;
% plot(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{motif_num}.FFvals_array_containsNaN(:,syl_pos));
% 
% % == raw data
% 
% rawdates=cell2mat(AllDays_RawDatStruct{1}.data.(syl)(:,
% plot(cell2mat(AllDays_RawDatStruct{1}.data.(syl)(:,1)))


%% PLOT VARIOUS EXPLORATORY THINGS ABOUT MOTIFS

NumMotifRendsWanted=2; % e.g. if 2 then will collect first 2 renditions of the motif and correlate them.

% === Are within song-values correlated?
% given a syllable in motif rend 1, how correlated is to the same syllable in
% motif 2? to other syls in motif 2?

% COLLECT ALL DATA - one line for motif rend 1, another for motif rend 2.
for i=1:NumMotifs;
    
    % === COLLECT all data from songs with at least two renditions of this
    % motif
    unique_songs=unique(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.Song_datenum_vals(:));
    c=1; % counter, across all songs with enough motifs

    for ii=1:length(unique_songs);
        songnum=unique_songs(ii);
        inds_this_song=find(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.Song_datenum_vals==songnum);
        
        if length(inds_this_song)<NumMotifRendsWanted;
            % then this song does not have enough rends of this motif
            continue
        end
        
        % == collect the first and second motifs.
        % get inds in order (of motif occurance)
        tmp=nanmean(AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.SylPos_in_song(inds_this_song,:),2);
        
        [~, inds]=sort(tmp);
       
        inds_this_song=inds_this_song(inds); % in order of motif occurance
        
        % ---- TAKE ONLY THE FIRST N motifs (e.g. 2)
        inds_this_song=inds_this_song(1:NumMotifRendsWanted);
        
        % === COLLECT MOTIF DATA
        dat=AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.FFvals_array_containsNaN(inds_this_song,:);
        
        % put into output structure
        AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.FFvals_MotifRendCorrs(:,:,c)=dat;
        
        c=c+1;
        
    end
end


        
% ==== PLOT
plot_cols=lt_make_plot_colors(NumMotifRendsWanted,0,0);

for i=1:NumMotifs;
    motif_syls=MotifList{i};
    
    for ii=1:length(motif_syls);
        
        syl=motif_syls{ii};
        
        lt_figure; title(syl); hold on;
        
        % === PLOT
        hplot=[];
        FFvals=[];
        for j=1:NumMotifRendsWanted;
            
            % === DAT
            FFvals{j}=AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{i}.FFvals_MotifRendCorrs(j,ii,:);
            FFvals{j}=squeeze(FFvals{j});
            
            hplot(j)=plot(FFvals{j}, 'Color',plot_cols{j});
            
        end
        
        % === PLOT SCATTER between 1st and 2nd motifs
        lt_figure; hold on; 
        
        lt_regress(FFvals{2}, FFvals{1}, 1);
        
        grid on; box on
        title(syl);
        xlabel('motif 1');
        ylabel('motif 2');
      
        Ylim=ylim;
        line([Ylim(1) Ylim(2)], [Ylim(1) Ylim(2)])

        legend(hplot,['Rendition num: ' [num2str(1:NumMotifRendsWanted)]])
    end
end




       







%% CORRELATIONS

motif_num=1;


FFvals=AllDays_PlotLearning.Correlations2.data_by_motifs.Baseline.motif_num{motif_num}.FFvals_array_containsNaN;


% Remove nans
rows_to_remove=sum(isnan(FFvals),2);

FFvals_nonan=FFvals;
FFvals_nonan(logical(rows_to_remove),:)=[];

rho=corr(FFvals_nonan);


lt_figure; hold on;
imagesc(rho);
colormap('gray')

rho(:,4)





%% PLOT RENDITIONS 



%% PARAMS

AllSylFields=fieldnames(AllDays_StructStatsStruct.IndivSyls);
BaselineDays=Params.SeqFilter.BaselineDays;


%% SUBTRACT DAY MEAN FOR EACH DAY

% all renditions, get z-score of feature vector relative to the
% mean/std feature vector of that day.
syllist=AllSylFields;
NumDays=datenum(Params.SeqFilter.LastDay)-datenum(Params.SeqFilter.FirstDay)+1;

for i=1:length(syllist);
    syl=syllist{i};
    
    % for each day, subtract mean of that day
    for ii=1:NumDays;
        
        inds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,ii);
        FVdayvals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(inds,:);
        FVdaymean=mean(FVdayvals);
        FVdaystd=std(FVdayvals);
        
        FVdayMinusMean=(FVdayvals-repmat(FVdaymean,size(FVdayvals,1),1));
        FVdayzscores=(FVdayvals-repmat(FVdaymean,size(FVdayvals,1),1))./repmat(FVdaystd,size(FVdayvals,1),1);
        
        
        AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(inds,:)=FVdayzscores;
        AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_MinusDayMean(inds,:)=FVdayMinusMean;
        
    end
end



%% PLOT ALL RENDITIONS TO VISUALIZE CORRELATIONS, ETC
for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder);

    syllist=Params.SeqFilter.SylLists.FieldsInOrder{j};

% 1) NOT SMOOTHED
figure; hold on;
title('All renditions - pitch vs. time');

PlotColors=lt_make_plot_colors(length(syllist),0);

for i=1:length(syllist);
    syl=syllist{i};
    
    
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(:,2);
   
      % PLOT
      % color based on relationship to targ syl.
      if  any(ismember(Params.SeqFilter.SylLists.TargetSyls,syl)); % target syl
              hplot(i)=plot(times.FinalValue,vals,'-k');
      elseif any(ismember(Params.SeqFilter.SylLists.SylsSame,syl)); % similar syl
              hplot(i)=plot(times.FinalValue,vals,'.','Color',rand*([0 1 0]));
      elseif any(ismember(Params.SeqFilter.SylLists.SylsDifferent,syl)); % similar syl
              hplot(i)=plot(times.FinalValue,vals,'.','Color',rand*([1 0 0]));
      end
              
end

legend(hplot,syllist)

% 2) SMOOTHED
smbin=5; % rends to smooth over
figure; hold on;
title(['Smoothed (bin: ' num2str(smbin) ' rends): Pitch vs. time for all syls']);

for i=1:length(syllist);
    syl=syllist{i};
    
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    times_sm=lt_running_stats(times.FinalValue,smbin);
    
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(:,2);
    vals_sm=lt_running_stats(vals,smbin);
   
      % PLOT
      % color based on relationship to targ syl.
      if  any(ismember(Params.SeqFilter.SylLists.TargetSyls,syl)); % target syl
              hplot(i)=plot(times_sm.Mean,vals_sm.Mean,'-ok');
      elseif any(ismember(Params.SeqFilter.SylLists.SylsSame,syl)); % similar syl
              hplot(i)=plot(times_sm.Mean,vals_sm.Mean,'-o','Color',rand*([0 1 0]));
      elseif any(ismember(Params.SeqFilter.SylLists.SylsDifferent,syl)); % similar syl
              hplot(i)=plot(times_sm.Mean,vals_sm.Mean,'-o','Color',rand*([1 0 0]));
      end

end

legend(hplot,syllist)


% 3) Z-SCORED, same as above
% smoothed
smbin=5; % rends to smooth over
figure; hold on;
title(['Smoothed (bin: ' num2str(smbin) ' rends): Pitch vs. time for all syls']);
hfig=[];
for i=1:length(syllist);
    hfig(i)=subplot(length(syllist),1,i); hold on;
    syl=syllist{i};
    
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    times_sm=lt_running_stats(times.FinalValue,smbin);
    
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(:,2);
    vals_sm=lt_running_stats(vals,smbin);
   
      % PLOT
      % color based on relationship to targ syl.
      if  any(ismember(Params.SeqFilter.SylLists.TargetSyls,syl)); % target syl
              hplot=plot(times_sm.Mean,vals_sm.Mean,'.k');
      elseif any(ismember(Params.SeqFilter.SylLists.SylsSame,syl)); % similar syl
              hplot=plot(times_sm.Mean,vals_sm.Mean,'.','Color',rand*([0 1 0]));
      elseif any(ismember(Params.SeqFilter.SylLists.SylsDifferent,syl)); % similar syl
              hplot=plot(times_sm.Mean,vals_sm.Mean,'.','Color',rand*([1 0 0]));
      end

legend(hplot,syllist{i});
end

linkaxes(hfig,'x');

end


%% Get Correlation matrix by having song as a single value - 
% take advatange of fact that all renditions of all syllable times in a
% given song have exactly the same time values

% ORGANIZE ALL DATA INTO A MATRIX - each column is a syl, each row is a
% timepoint, and each element is average value in that timepoint (=song)

% convert to inds
if ischar(DaysWanted);
    if DaysWanted=='baseline';
        DayIndsWanted=BaselineDays;
    else
        disp('Problem, not sure what days wanted for corr matrix');
    end
elseif isnumeric(DaysWanted);
    DayIndsWanted=DaysWanted;
end



% get all times - these times required to contain all syls
syl=Params.SeqFilter.SylLists.TargetSyls{1}; % can use any syl to get times

% Only look at baseline days
bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DayIndsWanted); % inds of rends that are from baseline
times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum(bInds); % times for all rends
times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times

SongTimes=unique(times.FinalValue); % find all unique times (i.e. songs), in baseline

% First concatenate syl names of all motifs (in order [motif 1 motif 2];
SylsAllMotifsOrdered={};
for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    SylsAllMotifsOrdered=[SylsAllMotifsOrdered Params.SeqFilter.SylLists.FieldsInOrder{j}];
end


% For each song time, get values for syls (1 mean value)
c=1;
syllist=SylsAllMotifsOrdered; % all syls from all motifs in order,

for i=1:length(SongTimes);
    curtime=SongTimes(i);
    
    for ii=1:length(syllist);
        syl=syllist{ii};
        
        % -- Extract times and values
        bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DayIndsWanted);
        times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum(bInds); % times for all rends
        times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
        times=times.FinalValue;
        
        vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(bInds,2);
        
        % -- for each song time, get mean value
        % find values for that song
        inds=ismember(times,curtime);
        
        % Save info of how many rends in each song
        RendsPerSong(i,ii)=sum(inds);
        
        % if there are renditions, then get mean
        if any(inds)==0; % i.e. no value for this song
            MATRIX(c,ii)=nan;
        else
            MATRIX(c,ii)=mean(vals(inds));
        end
        
    end
    
    % check if any syls lacked data. if so, remove this row (by not
    % letting "c" add one more
    if any(isnan(MATRIX(c,:)));
        MATRIX(c,:)=[];
        disp(['THREW OUT SONG num (approx) ' num2str(c) ' because lacked data for at least one syl']);
        continue
    end
    
    c=c+1; % use this becuase some songs will throw out if does not have dat for all syls.
end

disp(['N (songs) used for corr: ' num2str(c-1)]);
disp(['N (songs) originally: ' num2str(length(SongTimes)) ]);


% GET CORRELATION MATRIX
[rho, p]=corr((MATRIX(:,:)));

% PLOT
figure; hold on;
colormap('jet');

subplot(2,1,1); hold on;
title(['Days ' num2str(DaysWanted) ': Correlation matrix (rho) betweens syls (mean val in song)']);
imagesc(rho,[-1 1]);
set(gca,'YTick',1:length(syllist));
set(gca,'YTickLabel',syllist);
set(gca,'XTick',1:length(syllist));
set(gca,'XTickLabel',syllist);

colorbar;

subplot(2,1,2); hold on;
title(['Days ' num2str(DaysWanted) ': matrix of p-values of correlation between syls (mean val in song)']);
imagesc(p,[-1 1]);
set(gca,'YTick',1:length(syllist));
set(gca,'YTickLabel',syllist);
set(gca,'XTick',1:length(syllist));
set(gca,'XTickLabel',syllist);

colorbar;




% TAKE 1D SLICE - correlations between target syl and others
for i=1:length(Params.SeqFilter.SylLists.TargetSyls); % for all targ syls
    targsyl=Params.SeqFilter.SylLists.TargetSyls{i};
    
    
    [rho, p]=corr(MATRIX); % gett correaltion matrix
    
    % find index of target syl.
    inds=strcmp(SylsAllMotifsOrdered,targsyl); % find where targ syl is in the corr matrix
    
    vals_rho=rho(:,inds); % extract column for targ syl vs. others
    vals_p=p(:,inds);
    
    figure; hold on;
    subplot(2,1,1); hold on;
    bar([vals_rho]);
    title(['Days ' num2str(DaysWanted) ' : rho of correlation between target (' targsyl ') and other syls']);
    set(gca,'XTick',1:length(SylsAllMotifsOrdered),'XTickLabel',SylsAllMotifsOrdered);
    
    subplot(2,1,2); hold on;
    bar([vals_p]);
    title(['Days ' num2str(DaysWanted) ' : p-value of corr between target (' targsyl ') and other syls']);
    set(gca,'XTick',1:length(SylsAllMotifsOrdered),'XTickLabel',SylsAllMotifsOrdered);
    
end

% -- CORRELATIONS PREDICT LEARNING?





%% CORRELATION FOR EACH DAY - look at evolution over time 
% OBSOLETE - USE ABOVE IN A FOR LOOP

%% GET CORRELATION MATRIX BY DIVIDING DATA USING TIME BINS
% will have a sliding time window - time windows in which both syls have
% enough samples will contribute to the correlation score.
HrBinSize=0.25; % size of bin, in hours
HrBinSlide=0.1; % amount to slide bin 
HrBinMinSamples=4; % minimum number of samples a bin must have (or it will be thrown out).







% Plot histogram of number of samples in all bins




