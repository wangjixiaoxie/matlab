%% LT 6/3/15 - Takes config file and desired notenum, and runs evtaf sim to get all detects, give hit infomration and freq

function [AllSongsData_toplot, AllData]=lt_check_hit_templ_freq_NoLabelRequired(Params)
% Instructions:
% Run in day folder. Only works for data collected with evtafv4 (since used
% notenum)


% Inputs, example:
% Params.batch='batch.rand.keep';
% Params.config='/bluejay4/lucas/birds/rd23gr89/060315_SeqDepPitch_durWN_day1/config.evconfig2'; 
% Params.NoteNum_to_plot=2; % for the note you want to analyze

% Outputs:
% AllSongsData_toplot - trial by trial for the target notenum only.

% TO DO: 
% does not do catch trial - will label all as hits, even if catch trial. 


%% -- Run evtafsim and save information of all songs
saveON=0;
AllData=lt_context_ExtractDayData(Params, saveON);


%%  TAKE ALLDATA AND extract trial by trial info

% == Get parameters:
% today's date from directory name
[~, ~, date, ~]=lt_get_birdname_date_from_dir(1);
firstday=date{2}; % i.e. today

NumSongs=length(AllData); % num songs

% ==== Collect all data from all songs to trial-by-trial, also figureing
% what what were actual hits vs. just detects.
AllSongsDataMatrix=[];

for i=1:NumSongs;
    
    n=size(AllData(i).data_OfflineAllDetects,1); % how many datapoint in this song.
    
    if n>0; % continue if this song has data.
        
        % -- 1) Look at online data to figure out which of those offline
        % renditions are hits vs. escapes.  Will determine which offline trials correspond
        % to online trials by comparing ttimes
        ttime_list_offline=AllData(i).data_OfflineAllDetects(:,1);
        
        online_inds_that_are_hit_LOGICALS=zeros(n,1); % prepare, all start off as not triggers. will modify to 1, those that are triggers
        online_inds_that_are_hit=[];
        
        if ~isempty(AllData(i).data_OnlineTrigs.ttimes);
            NumTrigs=size(AllData(i).data_OnlineTrigs.ttimes,1);
            
            for j=1:NumTrigs; % for all datapouints
                ttime_online=AllData(i).data_OnlineTrigs.ttimes(j);
                trignote=AllData(i).data_OnlineTrigs.trignotes(j);
                
                % -- which offline ttimes are closes to this online ttime? (could be multiple offline ttimes);
                
                % Note: empirically seems as if offline ttimes are never more
                % than 5ms diff from online trig times. Given that multiple
                % trigs never really happen within trig refract of <10ms, will
                % use 10ms as the value below which I say offline and online
                % ttimes are the same.  Above this value there are sometimes
                % very different values, indicating that there was an online
                % hit that was not detected online.  I throw these out. Low
                % fraction (2/289; 7/444, on 2 days for rd23gr89). Seems to be related to WN beign on that
                % sylalble?
                
                % If using template offline that is different from template
                % used online, then seems to see almost all time
                % differences be <25. So then would use 25 instead of 10.
                
                online_inds_that_are_hit=[online_inds_that_are_hit; find(abs(ttime_list_offline-ttime_online)<10)];
            end
        end
        
        online_inds_that_are_hit=unique(online_inds_that_are_hit);
        online_inds_that_are_hit_LOGICALS(online_inds_that_are_hit)=1;
        
        
        % -- 2) absolute time of song
        tmp=lt_convert_EventTimes_to_RelTimes(firstday,AllData(i).datenum_filestart_SecResol);
        %         [~, tmp]=lt_convert_datenum_to_hour(AllData(i).datenum_filestart_SecResol);
        timeofsong_days=tmp.FinalValue;
        
        % --- 3) Catch status of each trial
        
        % -- 4) Collect offline data
        AllSongsDataMatrix=[AllSongsDataMatrix; i*ones(n,1), timeofsong_days*ones(n,1), AllData(i).data_OfflineAllDetects, online_inds_that_are_hit_LOGICALS];
    end
end


% ====== Make sure data are in chronological order
[~, inds]=sort(AllSongsDataMatrix(:,2)); % list of times

AllSongsDataMatrix=AllSongsDataMatrix(inds,:);


%% Sort out just datapoints that are detects of desired note
AllSongsData_toplot=AllSongsDataMatrix(AllSongsDataMatrix(:,4)==Params.NoteNum_to_plot,:);


%% ==== PLOT DATA
RunBin=20;

% ====== 1) ALL RENDS (by time)
lt_figure; hold on;
title(['FF of all detects of notenum: ' num2str(Params.NoteNum_to_plot) ' (by time)']);
xlabel('time (hours)');
ylabel('FF (hz)');

% all rends
X=AllSongsData_toplot(:,2); % time values
% convert times from days to hours
X=(X-1)*24;
Y=AllSongsData_toplot(:,5); % freq values
plot(X, Y,'o','MarkerFaceColor','b');

% get running averages
X_sm=lt_running_stats(X,RunBin);
Y_sm=lt_running_stats(Y,RunBin);

% overlay online hits
IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
X=AllSongsData_toplot(IndsOfTriggeredRends,2); % only those triggered
% convert times from days to hours
X=(X-1)*24;

Y=AllSongsData_toplot(IndsOfTriggeredRends,5);
plot(X,Y,'or','MarkerFaceColor','r');

% overlay running avg
shadedErrorBar(X_sm.Mean,Y_sm.Mean,Y_sm.SEM,{'-k'},1);


% ===== 2) RUNNING PERCENTILES
X=AllSongsData_toplot(:,2); % time values
% convert times from days to hours
X=(X-1)*24;
Y=AllSongsData_toplot(:,5); % freq values

Y_running=lt_running_stats(Y,RunBin);
X_running=lt_running_stats(X,RunBin);
% -- PLOT 30 and 70 tiles
lt_figure; hold on;
title(['30th and 70th percentiles, binsize: ' num2str(RunBin) ' rends']);
xlabel('median time of bin (hours)');

% 30th
lt_plot(X_running.Median, Y_running.prctiles_5_30_50_70_95(:,2),{'Color', 'b'});
lt_plot(X_running.Median, Y_running.prctiles_5_30_50_70_95(:,4), {'Color','r'});


% ======= 3) Running average of hit rate
% 1) time points of all data
X=AllSongsData_toplot(:,2); % time values
% convert times from days to hours
X=(X-1)*24;

% 2) timepoints of only hit data
IndsOfTriggeredRends=find(AllSongsData_toplot(:,7)==1);
X_hit=AllSongsData_toplot(IndsOfTriggeredRends,2); % only those triggered
% convert times from days to hours
X_hit=(X_hit-1)*24;

% create bin edges using first and last timepoint (for histogram)
hist_edges=linspace(X(1), X(end), 10);

[X_hist]=histc(X, hist_edges);
[X_hit_hist]=histc(X_hit, hist_edges);

Running_Hit_Rate=X_hit_hist./X_hist;

lt_figure; hold on;
title('Running hit rate (catch songs count as no hits)');
xlabel('timepoint');
ylabel('hit rate');

bar(hist_edges, Running_Hit_Rate, 'histc');






%% SAVE
try 
    cd lt_check_hit_templ_freq_NoLabelRequired
catch err
    mkdir lt_check_hit_templ_freq_NoLabelRequired
    cd lt_check_hit_templ_freq_NoLabelRequired
end

tstamp=lt_get_timestamp(0);

mkdir(tstamp);

cd(tstamp)

save('AllSongsData_toplot','AllSongsData_toplot');
save('AllData','AllData')

lt_save_all_figs;

cd ../
cd ..



