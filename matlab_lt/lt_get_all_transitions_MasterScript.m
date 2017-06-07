%% 6/2/17 - LT, written quickly for Chrstian, for quick analysis of sequence learning.

%% first extract for each day
batchf = 'batch.labeled.all';
nameofbird = 'pu11wh87';
folder_phrase = 'SeqDepPitchShift2';
datestring = '21Nov2014';
phrase_filename = '';


lt_get_all_transition_probabilities_FUNCTION(batchf, nameofbird, ...
    folder_phrase,phrase_filename,datestring);

%% compile across days

% ================= PARAMS
clear all; close all;
SaveFiles = {'pu11wh87_21Nov2014_.mat', ...
    'pu11wh87_22Nov2014_.mat'};

TransitionsToPlot = {'ba', 'bb'};

FirstDay = '21Nov2014';
BinSizeSongs = 5;

plotRendsPerSong = 0; % if 0, then prob of transition per song

% ================== RUNS
lt_figure; hold on;
plotcols = lt_make_plot_colors(length(TransitionsToPlot), 0, 0);

for i=1:length(SaveFiles)
    
    tmp = load(SaveFiles{i});
  
    RendsPerSongAll = []; % trans x song;
    TotalRends = [];
    NumSongs = [];
    TimeOfSongAll = [];
    for ii=1:length(TransitionsToPlot)
       
        syl1 = TransitionsToPlot{ii}(1);
        syl2 = TransitionsToPlot{ii}(2);
        
        
        RendsPerSong = cell2mat(tmp.all_trans.(syl1).transition_to_.(syl2).number_of_instances_persong);
        TimeOfSong = [];
        
        for j=1:length(tmp.all_trans.syl_times)
            
            fntime = median(tmp.all_trans.syl_times{j});
            fntime = lt_convert_EventTimes_to_RelTimes(FirstDay, fntime);
            
            TimeOfSong = [TimeOfSong fntime.FinalValue];
        
        end
        
        if plotRendsPerSong ==1
        plot(TimeOfSong, RendsPerSong, '.', 'Color', plotcols{ii});
        
        % --- running average
        rendsRunning = lt_running_stats(RendsPerSong, BinSizeSongs);
        timeRunning = lt_running_stats(TimeOfSong, BinSizeSongs);
        
        shadedErrorBar(timeRunning.Mean, rendsRunning.Mean, rendsRunning.SEM, {'Color', plotcols{ii}}, 1);
        end
        
        % -- get total rends
        TotalRends = [TotalRends sum(RendsPerSong)];
        NumSongs = [NumSongs length(RendsPerSong)];
        RendsPerSongAll = [RendsPerSongAll; RendsPerSong];
        TimeOfSongAll = [TimeOfSongAll; TimeOfSong];
    end
    assert(all(diff(NumSongs)==0), 'asfasd');
    
    if plotRendsPerSong ==0
       for ii=1:length(TransitionsToPlot)
          ProbPerSong = RendsPerSongAll(ii, :)./sum(RendsPerSongAll, 1); 
          
          plot(TimeOfSongAll(ii, :), ProbPerSong, '-o', 'Color', plotcols{ii});
           
       end
    end
       
    
    % --- plot total rends
    for ii=1:length(TotalRends)
       X =  floor(TimeOfSong)+0.9 + 0.05*ii;
       if plotRendsPerSong==1
       Y = TotalRends(ii)/NumSongs(1); % mean number of rends per song
       elseif plotRendsPerSong==0
       Y = TotalRends(ii)/sum(TotalRends); % mean number of rends per song
       end
       lt_plot(X, Y, {'Color', plotcols{ii}});
    end

end

